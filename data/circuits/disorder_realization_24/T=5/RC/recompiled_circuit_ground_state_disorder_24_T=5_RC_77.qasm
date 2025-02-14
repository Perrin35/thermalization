OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.71819031) q[0];
sx q[0];
rz(-3.1388404) q[0];
sx q[0];
rz(1.0116853) q[0];
rz(-2.6818795) q[1];
sx q[1];
rz(-1.3016394) q[1];
sx q[1];
rz(0.17170061) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3637488) q[0];
sx q[0];
rz(-0.96723377) q[0];
sx q[0];
rz(1.6700755) q[0];
rz(-pi) q[1];
rz(2.6284211) q[2];
sx q[2];
rz(-2.9193239) q[2];
sx q[2];
rz(-1.5813511) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57331177) q[1];
sx q[1];
rz(-2.3900095) q[1];
sx q[1];
rz(-2.4922396) q[1];
rz(-pi) q[2];
rz(-3.0875077) q[3];
sx q[3];
rz(-1.4776609) q[3];
sx q[3];
rz(-2.2786811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1823938) q[2];
sx q[2];
rz(-2.6540519) q[2];
sx q[2];
rz(-1.7215151) q[2];
rz(-1.9567664) q[3];
sx q[3];
rz(-1.607837) q[3];
sx q[3];
rz(2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(3.0576393) q[0];
sx q[0];
rz(-1.5204484) q[0];
sx q[0];
rz(0.49348304) q[0];
rz(2.3656942) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(-0.50481558) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1349943) q[0];
sx q[0];
rz(-2.2033443) q[0];
sx q[0];
rz(-2.452144) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4846862) q[2];
sx q[2];
rz(-2.5599766) q[2];
sx q[2];
rz(1.5145375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9318092) q[1];
sx q[1];
rz(-1.1714913) q[1];
sx q[1];
rz(-0.65815355) q[1];
x q[2];
rz(-0.22498954) q[3];
sx q[3];
rz(-0.27249042) q[3];
sx q[3];
rz(-0.30411965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33727553) q[2];
sx q[2];
rz(-1.3108459) q[2];
sx q[2];
rz(-2.6709225) q[2];
rz(-0.63052952) q[3];
sx q[3];
rz(-2.1157406) q[3];
sx q[3];
rz(1.7414909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0640963) q[0];
sx q[0];
rz(-3.016576) q[0];
sx q[0];
rz(-2.4858544) q[0];
rz(-2.8302622) q[1];
sx q[1];
rz(-0.89404023) q[1];
sx q[1];
rz(-3.0701367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27256672) q[0];
sx q[0];
rz(-1.9211968) q[0];
sx q[0];
rz(-0.053215543) q[0];
rz(-pi) q[1];
rz(-2.803903) q[2];
sx q[2];
rz(-1.2215541) q[2];
sx q[2];
rz(1.8943292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25117043) q[1];
sx q[1];
rz(-0.4931207) q[1];
sx q[1];
rz(-1.7095196) q[1];
x q[2];
rz(1.7830332) q[3];
sx q[3];
rz(-2.0925412) q[3];
sx q[3];
rz(2.9468342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91325703) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(0.060001686) q[2];
rz(1.2126728) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(2.3771299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5948828) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(-0.57719624) q[0];
rz(2.3120841) q[1];
sx q[1];
rz(-2.0894158) q[1];
sx q[1];
rz(2.3037691) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71815495) q[0];
sx q[0];
rz(-0.95036794) q[0];
sx q[0];
rz(-1.6654832) q[0];
rz(1.6610386) q[2];
sx q[2];
rz(-2.549516) q[2];
sx q[2];
rz(-0.83276487) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41552222) q[1];
sx q[1];
rz(-0.94968098) q[1];
sx q[1];
rz(-2.6776621) q[1];
rz(-pi) q[2];
rz(2.1076848) q[3];
sx q[3];
rz(-0.75124012) q[3];
sx q[3];
rz(0.43207016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.025658) q[2];
sx q[2];
rz(-1.6178774) q[2];
sx q[2];
rz(-1.9177829) q[2];
rz(0.4661679) q[3];
sx q[3];
rz(-0.40188447) q[3];
sx q[3];
rz(0.60552067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25662988) q[0];
sx q[0];
rz(-1.7054568) q[0];
sx q[0];
rz(3.0404941) q[0];
rz(2.7283607) q[1];
sx q[1];
rz(-0.44094545) q[1];
sx q[1];
rz(-2.0535927) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1285217) q[0];
sx q[0];
rz(-0.59935843) q[0];
sx q[0];
rz(-2.1764663) q[0];
x q[1];
rz(-1.1674983) q[2];
sx q[2];
rz(-2.3638862) q[2];
sx q[2];
rz(0.68439999) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1428296) q[1];
sx q[1];
rz(-1.3266449) q[1];
sx q[1];
rz(1.9119838) q[1];
x q[2];
rz(-3.0524714) q[3];
sx q[3];
rz(-2.4702419) q[3];
sx q[3];
rz(-1.7969839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92855144) q[2];
sx q[2];
rz(-0.88625675) q[2];
sx q[2];
rz(1.586033) q[2];
rz(-2.0689615) q[3];
sx q[3];
rz(-1.6173247) q[3];
sx q[3];
rz(-3.0090289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17305408) q[0];
sx q[0];
rz(-2.2569077) q[0];
sx q[0];
rz(-1.7065077) q[0];
rz(0.75812078) q[1];
sx q[1];
rz(-2.4393612) q[1];
sx q[1];
rz(3.039956) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1124737) q[0];
sx q[0];
rz(-0.85364193) q[0];
sx q[0];
rz(-0.6458592) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2287057) q[2];
sx q[2];
rz(-2.5882571) q[2];
sx q[2];
rz(-1.1007512) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96856252) q[1];
sx q[1];
rz(-0.30245879) q[1];
sx q[1];
rz(-0.36424251) q[1];
rz(-pi) q[2];
rz(-2.6061344) q[3];
sx q[3];
rz(-1.0046008) q[3];
sx q[3];
rz(-0.2406075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8018084) q[2];
sx q[2];
rz(-2.0399751) q[2];
sx q[2];
rz(1.8088809) q[2];
rz(2.0594275) q[3];
sx q[3];
rz(-2.3887631) q[3];
sx q[3];
rz(1.4235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4448755) q[0];
sx q[0];
rz(-1.8901905) q[0];
sx q[0];
rz(2.3984997) q[0];
rz(-0.7630868) q[1];
sx q[1];
rz(-0.63947314) q[1];
sx q[1];
rz(2.2230164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6800623) q[0];
sx q[0];
rz(-1.5627075) q[0];
sx q[0];
rz(-1.5458406) q[0];
rz(2.0430365) q[2];
sx q[2];
rz(-1.8794606) q[2];
sx q[2];
rz(-0.93738467) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1700701) q[1];
sx q[1];
rz(-0.60107175) q[1];
sx q[1];
rz(0.55723377) q[1];
rz(-1.6878456) q[3];
sx q[3];
rz(-2.2915386) q[3];
sx q[3];
rz(-0.085214867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4172198) q[2];
sx q[2];
rz(-1.1823187) q[2];
sx q[2];
rz(-0.43583885) q[2];
rz(-3.0271652) q[3];
sx q[3];
rz(-2.8947713) q[3];
sx q[3];
rz(2.54336) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33335394) q[0];
sx q[0];
rz(-2.5840608) q[0];
sx q[0];
rz(2.5575141) q[0];
rz(-2.7059879) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(1.6216507) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0348957) q[0];
sx q[0];
rz(-0.55904065) q[0];
sx q[0];
rz(2.1513272) q[0];
rz(0.21824117) q[2];
sx q[2];
rz(-1.9710961) q[2];
sx q[2];
rz(-1.2117653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78306373) q[1];
sx q[1];
rz(-1.7768246) q[1];
sx q[1];
rz(1.5611783) q[1];
x q[2];
rz(-2.5831843) q[3];
sx q[3];
rz(-0.30277751) q[3];
sx q[3];
rz(1.4608135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0674151) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(1.1161067) q[2];
rz(1.762278) q[3];
sx q[3];
rz(-2.0383056) q[3];
sx q[3];
rz(3.0232159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.24935687) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(0.78999162) q[0];
rz(3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(2.7008609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249464) q[0];
sx q[0];
rz(-1.2727203) q[0];
sx q[0];
rz(-1.9008725) q[0];
x q[1];
rz(-1.1665383) q[2];
sx q[2];
rz(-2.031293) q[2];
sx q[2];
rz(-1.3590036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3778119) q[1];
sx q[1];
rz(-1.0965938) q[1];
sx q[1];
rz(0.18486102) q[1];
x q[2];
rz(-1.3814028) q[3];
sx q[3];
rz(-1.5301068) q[3];
sx q[3];
rz(0.19408801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8792087) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(-0.66514307) q[2];
rz(2.3696259) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36122286) q[0];
sx q[0];
rz(-0.64798111) q[0];
sx q[0];
rz(1.004647) q[0];
rz(-1.7338344) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
sx q[1];
rz(-1.8515324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125975) q[0];
sx q[0];
rz(-0.83461232) q[0];
sx q[0];
rz(-0.33552977) q[0];
rz(-2.9384841) q[2];
sx q[2];
rz(-2.1611193) q[2];
sx q[2];
rz(-0.029435722) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0755578) q[1];
sx q[1];
rz(-1.7820638) q[1];
sx q[1];
rz(2.8610364) q[1];
x q[2];
rz(2.9674888) q[3];
sx q[3];
rz(-2.1645567) q[3];
sx q[3];
rz(1.6792149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1349692) q[2];
sx q[2];
rz(-1.9776191) q[2];
sx q[2];
rz(0.25035614) q[2];
rz(-0.97744673) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(-1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95494315) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(0.74465887) q[1];
sx q[1];
rz(-0.79569334) q[1];
sx q[1];
rz(-2.4462499) q[1];
rz(-1.5581239) q[2];
sx q[2];
rz(-1.7141533) q[2];
sx q[2];
rz(-1.2319596) q[2];
rz(1.0741735) q[3];
sx q[3];
rz(-1.1840829) q[3];
sx q[3];
rz(0.89331762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
