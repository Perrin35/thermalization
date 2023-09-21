OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5157226) q[0];
sx q[0];
rz(-0.54870257) q[0];
sx q[0];
rz(-0.8843511) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(-1.5024827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2384773) q[0];
sx q[0];
rz(-1.7503386) q[0];
sx q[0];
rz(1.205501) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56514481) q[2];
sx q[2];
rz(-1.5663212) q[2];
sx q[2];
rz(-2.7067513) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2734387) q[1];
sx q[1];
rz(-1.8167129) q[1];
sx q[1];
rz(1.5032561) q[1];
rz(-2.2803335) q[3];
sx q[3];
rz(-1.2951295) q[3];
sx q[3];
rz(2.0303126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(0.65594977) q[2];
rz(-1.2077228) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(-0.22878376) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-0.0072335009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4097071) q[0];
sx q[0];
rz(-1.1792372) q[0];
sx q[0];
rz(-3.0407228) q[0];
x q[1];
rz(-2.8480808) q[2];
sx q[2];
rz(-1.268317) q[2];
sx q[2];
rz(0.40700618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.05939535) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(2.5897964) q[1];
x q[2];
rz(-2.238027) q[3];
sx q[3];
rz(-1.9575319) q[3];
sx q[3];
rz(-2.8500593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(3.0200322) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4678629) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(-1.3695705) q[0];
rz(-1.9000152) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(-2.8799768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152552) q[0];
sx q[0];
rz(-1.5528423) q[0];
sx q[0];
rz(-3.1319588) q[0];
x q[1];
rz(0.12446603) q[2];
sx q[2];
rz(-2.3111768) q[2];
sx q[2];
rz(0.24701961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6690327) q[1];
sx q[1];
rz(-1.4497888) q[1];
sx q[1];
rz(2.3585412) q[1];
rz(3.0858585) q[3];
sx q[3];
rz(-2.627943) q[3];
sx q[3];
rz(1.0518215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(2.938081) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(-0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(1.571636) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(1.2483695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49444775) q[0];
sx q[0];
rz(-1.5054504) q[0];
sx q[0];
rz(1.6676184) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20935697) q[2];
sx q[2];
rz(-1.7484089) q[2];
sx q[2];
rz(2.3043485) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9256546) q[1];
sx q[1];
rz(-0.92080322) q[1];
sx q[1];
rz(-1.4309806) q[1];
rz(-pi) q[2];
rz(-0.72327153) q[3];
sx q[3];
rz(-1.372882) q[3];
sx q[3];
rz(-2.7151782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1127597) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(-0.83703414) q[2];
rz(1.933243) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(-2.3560431) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(0.77520448) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(2.2391589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8474903) q[0];
sx q[0];
rz(-0.93203629) q[0];
sx q[0];
rz(-1.5139447) q[0];
rz(-1.3257741) q[2];
sx q[2];
rz(-1.1282215) q[2];
sx q[2];
rz(-2.4545836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9296226) q[1];
sx q[1];
rz(-2.3096482) q[1];
sx q[1];
rz(-2.3421939) q[1];
x q[2];
rz(1.5794472) q[3];
sx q[3];
rz(-0.31673613) q[3];
sx q[3];
rz(2.1474311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9465785) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(1.9449332) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.8701575) q[3];
sx q[3];
rz(1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(-2.6384171) q[0];
rz(-1.4563837) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-2.9398289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87529463) q[0];
sx q[0];
rz(-1.6365956) q[0];
sx q[0];
rz(1.468303) q[0];
x q[1];
rz(2.6642338) q[2];
sx q[2];
rz(-1.1672033) q[2];
sx q[2];
rz(1.3178283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7355431) q[1];
sx q[1];
rz(-0.65945259) q[1];
sx q[1];
rz(0.72592782) q[1];
rz(-pi) q[2];
rz(-1.1062578) q[3];
sx q[3];
rz(-2.2482276) q[3];
sx q[3];
rz(0.57074742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(2.8743437) q[2];
rz(-0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(-0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8478407) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-3.0840432) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4276886) q[0];
sx q[0];
rz(-1.9644992) q[0];
sx q[0];
rz(-2.4967771) q[0];
rz(0.96037453) q[2];
sx q[2];
rz(-0.27563169) q[2];
sx q[2];
rz(-0.20197091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6616933) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(-2.2366621) q[1];
rz(-pi) q[2];
rz(-1.1293291) q[3];
sx q[3];
rz(-0.55674508) q[3];
sx q[3];
rz(2.7260775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0223579) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(0.38267246) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(-0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.7169749) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(2.8714645) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547202) q[0];
sx q[0];
rz(-0.76612872) q[0];
sx q[0];
rz(0.92673577) q[0];
rz(-pi) q[1];
rz(0.59992744) q[2];
sx q[2];
rz(-1.6919961) q[2];
sx q[2];
rz(-1.3231414) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24208454) q[1];
sx q[1];
rz(-0.71166066) q[1];
sx q[1];
rz(0.74739026) q[1];
x q[2];
rz(1.6007401) q[3];
sx q[3];
rz(-0.91120126) q[3];
sx q[3];
rz(-2.3823882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55390629) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(-1.2109057) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07638409) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(2.912345) q[0];
rz(-0.30300888) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(1.4607666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0871353) q[0];
sx q[0];
rz(-1.5532877) q[0];
sx q[0];
rz(0.28355916) q[0];
x q[1];
rz(-0.42713366) q[2];
sx q[2];
rz(-0.41296994) q[2];
sx q[2];
rz(0.6461179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2673805) q[1];
sx q[1];
rz(-1.0460209) q[1];
sx q[1];
rz(-1.5449779) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78482307) q[3];
sx q[3];
rz(-1.8514957) q[3];
sx q[3];
rz(2.7626038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.71172697) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(-0.7152043) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-0.33690548) q[0];
rz(0.20740549) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-2.7609603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.80523) q[0];
sx q[0];
rz(-0.21348937) q[0];
sx q[0];
rz(-1.2345033) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1516018) q[2];
sx q[2];
rz(-0.91772807) q[2];
sx q[2];
rz(2.4534015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8457348) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(2.9198398) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57772824) q[3];
sx q[3];
rz(-1.2168222) q[3];
sx q[3];
rz(2.3059394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(-1.1364737) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.3142746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3363591) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(-1.0271172) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(-1.3409875) q[2];
sx q[2];
rz(-1.8125712) q[2];
sx q[2];
rz(-0.3704091) q[2];
rz(2.9410578) q[3];
sx q[3];
rz(-1.1391098) q[3];
sx q[3];
rz(-1.1548635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];