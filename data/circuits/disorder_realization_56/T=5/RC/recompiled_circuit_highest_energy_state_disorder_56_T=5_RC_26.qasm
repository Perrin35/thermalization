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
rz(1.4752969) q[0];
sx q[0];
rz(1.8721606) q[0];
sx q[0];
rz(8.7526487) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(-2.3840005) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8920994) q[0];
sx q[0];
rz(-1.4612911) q[0];
sx q[0];
rz(3.0654869) q[0];
rz(2.5301928) q[2];
sx q[2];
rz(-0.58908236) q[2];
sx q[2];
rz(-3.0036199) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.014908286) q[1];
sx q[1];
rz(-1.1619301) q[1];
sx q[1];
rz(-1.3247299) q[1];
rz(2.4663062) q[3];
sx q[3];
rz(-2.3216341) q[3];
sx q[3];
rz(0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0932833) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(-1.7993571) q[2];
rz(1.7736769) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(1.5889408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9369478) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(2.192705) q[0];
rz(-3.0401547) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(0.92996517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0895871) q[0];
sx q[0];
rz(-1.4692592) q[0];
sx q[0];
rz(-2.9278281) q[0];
rz(0.50518121) q[2];
sx q[2];
rz(-1.564216) q[2];
sx q[2];
rz(-0.40019658) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3392433) q[1];
sx q[1];
rz(-1.7319458) q[1];
sx q[1];
rz(-1.4646962) q[1];
rz(-pi) q[2];
rz(-0.92949246) q[3];
sx q[3];
rz(-0.87021032) q[3];
sx q[3];
rz(1.6456347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0019504) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(-2.3993717) q[2];
rz(0.46418515) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6716229) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(-1.4416913) q[0];
rz(-0.16547671) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(-2.2742719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21737032) q[0];
sx q[0];
rz(-1.570756) q[0];
sx q[0];
rz(1.5713619) q[0];
rz(-pi) q[1];
rz(-0.27133743) q[2];
sx q[2];
rz(-2.0607161) q[2];
sx q[2];
rz(2.4957531) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9359303) q[1];
sx q[1];
rz(-1.8962269) q[1];
sx q[1];
rz(1.8673351) q[1];
rz(-pi) q[2];
rz(-2.5977449) q[3];
sx q[3];
rz(-2.2529054) q[3];
sx q[3];
rz(2.2206375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9637588) q[2];
sx q[2];
rz(-1.1744262) q[2];
sx q[2];
rz(-2.3826694) q[2];
rz(0.80398503) q[3];
sx q[3];
rz(-2.1920429) q[3];
sx q[3];
rz(2.4132531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8300962) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-0.48402825) q[0];
rz(-1.5361891) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(1.4253634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.037402) q[0];
sx q[0];
rz(-2.0811102) q[0];
sx q[0];
rz(-1.9602106) q[0];
x q[1];
rz(0.31230782) q[2];
sx q[2];
rz(-1.7877525) q[2];
sx q[2];
rz(1.7157451) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4933136) q[1];
sx q[1];
rz(-1.7910622) q[1];
sx q[1];
rz(-0.92695285) q[1];
x q[2];
rz(1.1953765) q[3];
sx q[3];
rz(-1.9945126) q[3];
sx q[3];
rz(-2.8590607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75590762) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(-0.29328406) q[2];
rz(-3.0934603) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(-0.3046681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83679477) q[0];
sx q[0];
rz(-1.7339107) q[0];
sx q[0];
rz(1.1037214) q[0];
rz(0.91148218) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(0.10890659) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2809362) q[0];
sx q[0];
rz(-2.0514279) q[0];
sx q[0];
rz(2.3970766) q[0];
rz(-pi) q[1];
rz(-1.3896732) q[2];
sx q[2];
rz(-2.0985773) q[2];
sx q[2];
rz(2.9865047) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2786633) q[1];
sx q[1];
rz(-0.68611523) q[1];
sx q[1];
rz(-1.8300232) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8725242) q[3];
sx q[3];
rz(-1.0142039) q[3];
sx q[3];
rz(0.8100141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64985046) q[2];
sx q[2];
rz(-1.7848585) q[2];
sx q[2];
rz(-2.8835127) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-0.93770599) q[3];
sx q[3];
rz(-0.55735731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7447164) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(-1.1599524) q[0];
rz(-0.57811919) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(0.32346183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1581381) q[0];
sx q[0];
rz(-2.0120512) q[0];
sx q[0];
rz(2.6652314) q[0];
rz(-0.49655621) q[2];
sx q[2];
rz(-1.7944031) q[2];
sx q[2];
rz(-1.8211435) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0486794) q[1];
sx q[1];
rz(-1.3447176) q[1];
sx q[1];
rz(-2.5552555) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46053912) q[3];
sx q[3];
rz(-1.2804739) q[3];
sx q[3];
rz(-3.0825305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9385927) q[2];
sx q[2];
rz(-0.77363571) q[2];
sx q[2];
rz(0.8482376) q[2];
rz(1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(-2.447824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(3.0271725) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(-2.2681336) q[0];
rz(-1.2365485) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(3.0580318) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3665584) q[0];
sx q[0];
rz(-1.8904443) q[0];
sx q[0];
rz(1.4421706) q[0];
rz(-pi) q[1];
rz(1.9651863) q[2];
sx q[2];
rz(-0.56496921) q[2];
sx q[2];
rz(2.625287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.04490964) q[1];
sx q[1];
rz(-2.2538141) q[1];
sx q[1];
rz(0.51814305) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7807981) q[3];
sx q[3];
rz(-1.8550145) q[3];
sx q[3];
rz(-2.7786215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7053335) q[2];
sx q[2];
rz(-1.5353563) q[2];
sx q[2];
rz(-1.7748888) q[2];
rz(-2.8691835) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(0.75850707) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2347539) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(0.93836623) q[0];
rz(2.3387108) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(-2.3209007) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9296449) q[0];
sx q[0];
rz(-1.4370455) q[0];
sx q[0];
rz(-1.9424214) q[0];
x q[1];
rz(0.16751473) q[2];
sx q[2];
rz(-1.9394839) q[2];
sx q[2];
rz(-2.9978254) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0487602) q[1];
sx q[1];
rz(-0.63221778) q[1];
sx q[1];
rz(-1.7507751) q[1];
x q[2];
rz(1.9441767) q[3];
sx q[3];
rz(-0.42694091) q[3];
sx q[3];
rz(-2.3641242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.20251033) q[2];
sx q[2];
rz(-0.15190092) q[2];
sx q[2];
rz(-1.067777) q[2];
rz(-1.1311401) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(-0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15602569) q[0];
sx q[0];
rz(-1.1556867) q[0];
sx q[0];
rz(0.40618968) q[0];
rz(-1.4093026) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(-1.0850151) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1472441) q[0];
sx q[0];
rz(-1.7450593) q[0];
sx q[0];
rz(-0.6263349) q[0];
rz(1.5662534) q[2];
sx q[2];
rz(-2.1811003) q[2];
sx q[2];
rz(1.4686327) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.93714) q[1];
sx q[1];
rz(-1.8984183) q[1];
sx q[1];
rz(-2.3069068) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75094838) q[3];
sx q[3];
rz(-0.72970684) q[3];
sx q[3];
rz(2.2580604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5558418) q[2];
sx q[2];
rz(-0.61183524) q[2];
sx q[2];
rz(-2.2108868) q[2];
rz(-1.7454923) q[3];
sx q[3];
rz(-1.3627108) q[3];
sx q[3];
rz(0.11955587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4569106) q[0];
sx q[0];
rz(-2.1645808) q[0];
sx q[0];
rz(1.8344301) q[0];
rz(2.7556509) q[1];
sx q[1];
rz(-1.7470876) q[1];
sx q[1];
rz(1.0345667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86901122) q[0];
sx q[0];
rz(-0.20244652) q[0];
sx q[0];
rz(-1.0342717) q[0];
x q[1];
rz(-1.3645646) q[2];
sx q[2];
rz(-2.0230556) q[2];
sx q[2];
rz(-2.3444676) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93232357) q[1];
sx q[1];
rz(-2.0210365) q[1];
sx q[1];
rz(2.2850014) q[1];
rz(-pi) q[2];
rz(2.3307269) q[3];
sx q[3];
rz(-1.4812143) q[3];
sx q[3];
rz(-0.60096622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4474386) q[2];
sx q[2];
rz(-2.4093781) q[2];
sx q[2];
rz(-2.7745957) q[2];
rz(-2.0224109) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(-1.5252349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3601396) q[0];
sx q[0];
rz(-1.1893138) q[0];
sx q[0];
rz(-1.0839533) q[0];
rz(-3.0837334) q[1];
sx q[1];
rz(-1.2670988) q[1];
sx q[1];
rz(1.7115464) q[1];
rz(-1.6690785) q[2];
sx q[2];
rz(-1.7429033) q[2];
sx q[2];
rz(1.8175816) q[2];
rz(0.69969768) q[3];
sx q[3];
rz(-1.8498004) q[3];
sx q[3];
rz(-0.39941272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
