OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(0.41419849) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.480455) q[0];
sx q[0];
rz(-1.9248795) q[0];
sx q[0];
rz(2.7829091) q[0];
x q[1];
rz(-2.1404999) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(2.9542838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84250433) q[1];
sx q[1];
rz(-0.43049225) q[1];
sx q[1];
rz(-1.4166142) q[1];
rz(2.8768086) q[3];
sx q[3];
rz(-1.6695108) q[3];
sx q[3];
rz(0.36987723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3027705) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(2.9717428) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-0.53952113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3300433) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(3.0898068) q[0];
rz(-pi) q[1];
rz(-2.4194854) q[2];
sx q[2];
rz(-1.8849843) q[2];
sx q[2];
rz(-2.9177641) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9425548) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(1.946279) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1971365) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(-0.32523793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(2.6932122) q[0];
rz(-1.386863) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019779531) q[0];
sx q[0];
rz(-0.67987961) q[0];
sx q[0];
rz(2.7133184) q[0];
x q[1];
rz(-1.0923549) q[2];
sx q[2];
rz(-1.1698327) q[2];
sx q[2];
rz(1.5329597) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2910034) q[1];
sx q[1];
rz(-0.91445078) q[1];
sx q[1];
rz(-2.7235051) q[1];
x q[2];
rz(-0.17685299) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(-2.0166486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(-0.57470542) q[2];
rz(-1.3556708) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280076) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(0.25948778) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703092) q[0];
sx q[0];
rz(-1.0571612) q[0];
sx q[0];
rz(-2.1156103) q[0];
x q[1];
rz(-0.71728431) q[2];
sx q[2];
rz(-1.6589763) q[2];
sx q[2];
rz(0.35536534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4983474) q[1];
sx q[1];
rz(-2.8162662) q[1];
sx q[1];
rz(-2.494032) q[1];
rz(-pi) q[2];
rz(1.4452403) q[3];
sx q[3];
rz(-1.351965) q[3];
sx q[3];
rz(0.20727508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(2.990492) q[2];
rz(-0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(-1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7933465) q[0];
sx q[0];
rz(-0.12013809) q[0];
sx q[0];
rz(-1.5915464) q[0];
rz(-0.30349489) q[2];
sx q[2];
rz(-1.3611088) q[2];
sx q[2];
rz(1.7393302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5835727) q[1];
sx q[1];
rz(-0.94545525) q[1];
sx q[1];
rz(-3.0261092) q[1];
rz(-1.2410774) q[3];
sx q[3];
rz(-1.1493249) q[3];
sx q[3];
rz(2.4333405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(2.8395555) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.1557895) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(-3.0715122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87175831) q[0];
sx q[0];
rz(-1.4176798) q[0];
sx q[0];
rz(1.7437177) q[0];
x q[1];
rz(-1.0191392) q[2];
sx q[2];
rz(-2.2746804) q[2];
sx q[2];
rz(0.64955074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1706714) q[1];
sx q[1];
rz(-2.1340003) q[1];
sx q[1];
rz(-1.7621653) q[1];
x q[2];
rz(2.0134986) q[3];
sx q[3];
rz(-3.0590995) q[3];
sx q[3];
rz(-0.024162956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.4403884) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(-2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550585) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(2.281718) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(0.0079356114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746349) q[0];
sx q[0];
rz(-1.8263706) q[0];
sx q[0];
rz(2.5711683) q[0];
x q[1];
rz(2.0769579) q[2];
sx q[2];
rz(-2.3828265) q[2];
sx q[2];
rz(-0.90492349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.36564) q[1];
sx q[1];
rz(-1.0675634) q[1];
sx q[1];
rz(-2.3535437) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.946297) q[3];
sx q[3];
rz(-2.7995031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(-2.725214) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(2.7767048) q[0];
rz(2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.6392802) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973328) q[0];
sx q[0];
rz(-1.5808006) q[0];
sx q[0];
rz(-0.25015932) q[0];
x q[1];
rz(-1.0566696) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(0.67509292) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8870526) q[1];
sx q[1];
rz(-2.2215543) q[1];
sx q[1];
rz(2.410143) q[1];
x q[2];
rz(3.0495166) q[3];
sx q[3];
rz(-1.3658938) q[3];
sx q[3];
rz(-0.74842194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-1.0650744) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(-0.43564963) q[0];
rz(-1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(-2.7246144) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57396736) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(0.91659878) q[0];
x q[1];
rz(-0.99879361) q[2];
sx q[2];
rz(-0.98888328) q[2];
sx q[2];
rz(1.7172161) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3560564) q[1];
sx q[1];
rz(-2.0691263) q[1];
sx q[1];
rz(2.9836125) q[1];
rz(-pi) q[2];
rz(3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(1.5537378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(0.22582516) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(-1.5243994) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(1.3226002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37227092) q[0];
sx q[0];
rz(-1.1936545) q[0];
sx q[0];
rz(-3.0110714) q[0];
x q[1];
rz(-2.0595423) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-0.79364712) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0364914) q[1];
sx q[1];
rz(-1.0579915) q[1];
sx q[1];
rz(0.75863104) q[1];
rz(-pi) q[2];
rz(2.2262276) q[3];
sx q[3];
rz(-1.4487106) q[3];
sx q[3];
rz(-1.9742427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(0.50921847) q[2];
sx q[2];
rz(-1.5794532) q[2];
sx q[2];
rz(3.0184359) q[2];
rz(1.3397459) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
