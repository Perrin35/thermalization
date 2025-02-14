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
rz(1.709197) q[0];
sx q[0];
rz(3.4442918) q[0];
sx q[0];
rz(11.533307) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(-1.6515825) q[1];
sx q[1];
rz(-2.9472247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5759566) q[0];
sx q[0];
rz(-2.4969184) q[0];
sx q[0];
rz(-2.8863532) q[0];
rz(-pi) q[1];
x q[1];
rz(1.871045) q[2];
sx q[2];
rz(-0.89648834) q[2];
sx q[2];
rz(0.44724321) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3857419) q[1];
sx q[1];
rz(-3.0839018) q[1];
sx q[1];
rz(2.2123446) q[1];
rz(-pi) q[2];
rz(2.5159149) q[3];
sx q[3];
rz(-1.2114269) q[3];
sx q[3];
rz(0.43454447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.73244652) q[2];
sx q[2];
rz(-1.1017841) q[2];
sx q[2];
rz(0.5644325) q[2];
rz(1.268528) q[3];
sx q[3];
rz(-0.0081491834) q[3];
sx q[3];
rz(-0.90659365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059435189) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(-0.91307688) q[0];
rz(-0.68436855) q[1];
sx q[1];
rz(-3.1412536) q[1];
sx q[1];
rz(2.2025462) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9371999) q[0];
sx q[0];
rz(-2.4449722) q[0];
sx q[0];
rz(-2.533769) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.026264391) q[2];
sx q[2];
rz(-2.677478) q[2];
sx q[2];
rz(0.26192203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0239183) q[1];
sx q[1];
rz(-0.54421762) q[1];
sx q[1];
rz(0.88440374) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7250118) q[3];
sx q[3];
rz(-0.41114488) q[3];
sx q[3];
rz(1.7257041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3188476) q[2];
sx q[2];
rz(-3.1343967) q[2];
sx q[2];
rz(2.2491271) q[2];
rz(1.5626296) q[3];
sx q[3];
rz(-0.024450863) q[3];
sx q[3];
rz(-1.310937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398294) q[0];
sx q[0];
rz(-2.6391397) q[0];
sx q[0];
rz(-0.40973642) q[0];
rz(-0.01092625) q[1];
sx q[1];
rz(-2.9210563) q[1];
sx q[1];
rz(1.3440557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4073915) q[0];
sx q[0];
rz(-2.0816861) q[0];
sx q[0];
rz(0.023500806) q[0];
rz(1.686269) q[2];
sx q[2];
rz(-1.6342548) q[2];
sx q[2];
rz(-2.7977914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8482159) q[1];
sx q[1];
rz(-2.8753706) q[1];
sx q[1];
rz(-1.8195527) q[1];
x q[2];
rz(0.009851708) q[3];
sx q[3];
rz(-1.6860699) q[3];
sx q[3];
rz(1.7719763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4695796) q[2];
sx q[2];
rz(-1.6894222) q[2];
sx q[2];
rz(-3.0984042) q[2];
rz(-2.0016661) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(-1.202701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4129055) q[0];
sx q[0];
rz(-0.40189704) q[0];
sx q[0];
rz(1.3380949) q[0];
rz(0.69522229) q[1];
sx q[1];
rz(-3.0137364) q[1];
sx q[1];
rz(-0.41740886) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25487568) q[0];
sx q[0];
rz(-2.4435197) q[0];
sx q[0];
rz(-1.3017734) q[0];
rz(-pi) q[1];
rz(2.874006) q[2];
sx q[2];
rz(-0.94174615) q[2];
sx q[2];
rz(-2.3305394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7898832) q[1];
sx q[1];
rz(-1.7725866) q[1];
sx q[1];
rz(-2.2191597) q[1];
rz(-pi) q[2];
rz(1.1379477) q[3];
sx q[3];
rz(-1.6995078) q[3];
sx q[3];
rz(-0.52279982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43856105) q[2];
sx q[2];
rz(-2.265354) q[2];
sx q[2];
rz(-1.38928) q[2];
rz(-2.321068) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(-2.4337721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97838068) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(2.1782844) q[0];
rz(2.9523201) q[1];
sx q[1];
rz(-0.015733868) q[1];
sx q[1];
rz(0.16545573) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236639) q[0];
sx q[0];
rz(-3.0907013) q[0];
sx q[0];
rz(2.5975157) q[0];
rz(2.8420429) q[2];
sx q[2];
rz(-0.89847696) q[2];
sx q[2];
rz(1.7393665) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.2442064) q[1];
sx q[1];
rz(-1.4546977) q[1];
sx q[1];
rz(2.3364324) q[1];
x q[2];
rz(-0.9521552) q[3];
sx q[3];
rz(-0.51636558) q[3];
sx q[3];
rz(-1.0960033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.307622) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(1.0597672) q[2];
rz(2.4506532) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(1.2714269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5237913) q[0];
sx q[0];
rz(-3.0484564) q[0];
sx q[0];
rz(1.5366489) q[0];
rz(0.45283428) q[1];
sx q[1];
rz(-3.1335399) q[1];
sx q[1];
rz(1.7013928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0217845) q[0];
sx q[0];
rz(-0.23384604) q[0];
sx q[0];
rz(2.3709557) q[0];
rz(-0.22162205) q[2];
sx q[2];
rz(-1.7852419) q[2];
sx q[2];
rz(0.947244) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3041087) q[1];
sx q[1];
rz(-2.9867801) q[1];
sx q[1];
rz(3.0938838) q[1];
x q[2];
rz(0.099248107) q[3];
sx q[3];
rz(-1.7407047) q[3];
sx q[3];
rz(2.3836294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77233934) q[2];
sx q[2];
rz(-0.99268308) q[2];
sx q[2];
rz(-3.0625647) q[2];
rz(-0.8655656) q[3];
sx q[3];
rz(-2.1871958) q[3];
sx q[3];
rz(-1.0237833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2257776) q[0];
sx q[0];
rz(-0.0053891698) q[0];
sx q[0];
rz(-0.22501568) q[0];
rz(0.30613884) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(-2.2711145) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3358094) q[0];
sx q[0];
rz(-0.085698232) q[0];
sx q[0];
rz(-0.58914574) q[0];
rz(-0.71484675) q[2];
sx q[2];
rz(-1.1230334) q[2];
sx q[2];
rz(3.062742) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6806769) q[1];
sx q[1];
rz(-0.80246882) q[1];
sx q[1];
rz(1.4867307) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87865307) q[3];
sx q[3];
rz(-2.5402398) q[3];
sx q[3];
rz(0.9200615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6656782) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(2.1062984) q[2];
rz(-0.78694844) q[3];
sx q[3];
rz(-0.042777177) q[3];
sx q[3];
rz(-0.27534819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.1110474) q[0];
sx q[0];
rz(-3.1304517) q[0];
sx q[0];
rz(0.025644843) q[0];
rz(-1.2643087) q[1];
sx q[1];
rz(-3.1176716) q[1];
sx q[1];
rz(0.69902507) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049533904) q[0];
sx q[0];
rz(-1.7932685) q[0];
sx q[0];
rz(0.27357863) q[0];
rz(-1.7412296) q[2];
sx q[2];
rz(-2.0336069) q[2];
sx q[2];
rz(-1.544726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1406724) q[1];
sx q[1];
rz(-0.39127054) q[1];
sx q[1];
rz(-2.9468582) q[1];
rz(-pi) q[2];
rz(0.91722971) q[3];
sx q[3];
rz(-1.1207657) q[3];
sx q[3];
rz(-0.91564028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25253025) q[2];
sx q[2];
rz(-0.74873304) q[2];
sx q[2];
rz(0.32386455) q[2];
rz(2.9845386) q[3];
sx q[3];
rz(-1.9040949) q[3];
sx q[3];
rz(0.15373716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5667628) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(0.59004849) q[0];
rz(0.7198965) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(2.2743026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7035422) q[0];
sx q[0];
rz(-0.64213404) q[0];
sx q[0];
rz(0.17583986) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9407104) q[2];
sx q[2];
rz(-2.5061766) q[2];
sx q[2];
rz(-0.70178343) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5979851) q[1];
sx q[1];
rz(-1.5026918) q[1];
sx q[1];
rz(-1.4582514) q[1];
rz(-pi) q[2];
rz(0.18901029) q[3];
sx q[3];
rz(-0.86652404) q[3];
sx q[3];
rz(-2.6420322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94816339) q[2];
sx q[2];
rz(-2.2488504) q[2];
sx q[2];
rz(-0.11290045) q[2];
rz(1.0490949) q[3];
sx q[3];
rz(-1.5821404) q[3];
sx q[3];
rz(-0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.076544) q[0];
sx q[0];
rz(-1.5600518) q[0];
sx q[0];
rz(1.6381868) q[0];
rz(-0.63649559) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(1.5719302) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51278567) q[0];
sx q[0];
rz(-0.11462002) q[0];
sx q[0];
rz(-0.8091457) q[0];
rz(-3.1384597) q[2];
sx q[2];
rz(-1.5701446) q[2];
sx q[2];
rz(-2.3359483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7074938) q[1];
sx q[1];
rz(-1.5705622) q[1];
sx q[1];
rz(-3.1386761) q[1];
rz(-pi) q[2];
rz(-1.4079553) q[3];
sx q[3];
rz(-2.0800839) q[3];
sx q[3];
rz(2.7247938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2703209) q[2];
sx q[2];
rz(-1.6297636) q[2];
sx q[2];
rz(0.16066571) q[2];
rz(-1.2284944) q[3];
sx q[3];
rz(-3.1077423) q[3];
sx q[3];
rz(0.20139774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(1.5184215) q[1];
sx q[1];
rz(-2.7802614) q[1];
sx q[1];
rz(-2.8526715) q[1];
rz(-1.2779909) q[2];
sx q[2];
rz(-1.59272) q[2];
sx q[2];
rz(-1.1772173) q[2];
rz(-2.0102228) q[3];
sx q[3];
rz(-1.7059587) q[3];
sx q[3];
rz(-0.95897924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
