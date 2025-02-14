OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1571932) q[0];
sx q[0];
rz(-2.0318883) q[0];
sx q[0];
rz(-1.00534) q[0];
rz(-2.4069064) q[1];
sx q[1];
rz(-0.35597304) q[1];
sx q[1];
rz(-2.2495143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1482502) q[0];
sx q[0];
rz(-2.3138232) q[0];
sx q[0];
rz(-1.4367075) q[0];
rz(-1.35688) q[2];
sx q[2];
rz(-2.5271673) q[2];
sx q[2];
rz(0.015903552) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6309384) q[1];
sx q[1];
rz(-2.0533685) q[1];
sx q[1];
rz(-1.4605182) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60601534) q[3];
sx q[3];
rz(-1.1258584) q[3];
sx q[3];
rz(-3.1387822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72508183) q[2];
sx q[2];
rz(-1.1434914) q[2];
sx q[2];
rz(-1.4408646) q[2];
rz(-2.1848988) q[3];
sx q[3];
rz(-1.1172833) q[3];
sx q[3];
rz(-0.24334894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956534) q[0];
sx q[0];
rz(-2.74701) q[0];
sx q[0];
rz(-1.12895) q[0];
rz(0.24540643) q[1];
sx q[1];
rz(-1.8281728) q[1];
sx q[1];
rz(0.66171563) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6889383) q[0];
sx q[0];
rz(-0.96315779) q[0];
sx q[0];
rz(1.2217962) q[0];
x q[1];
rz(-2.6241669) q[2];
sx q[2];
rz(-2.708507) q[2];
sx q[2];
rz(2.9703857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.052110869) q[1];
sx q[1];
rz(-0.27325892) q[1];
sx q[1];
rz(-3.0237126) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8850967) q[3];
sx q[3];
rz(-0.30499015) q[3];
sx q[3];
rz(1.5105607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5976065) q[2];
sx q[2];
rz(-0.49763766) q[2];
sx q[2];
rz(1.0428766) q[2];
rz(-1.4526224) q[3];
sx q[3];
rz(-2.2904604) q[3];
sx q[3];
rz(-2.0378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747028) q[0];
sx q[0];
rz(-0.68793982) q[0];
sx q[0];
rz(1.3091298) q[0];
rz(-0.4492999) q[1];
sx q[1];
rz(-0.6895014) q[1];
sx q[1];
rz(-1.7151394) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9096722) q[0];
sx q[0];
rz(-1.8350775) q[0];
sx q[0];
rz(2.8876165) q[0];
x q[1];
rz(1.6108247) q[2];
sx q[2];
rz(-0.42280254) q[2];
sx q[2];
rz(1.947248) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92628463) q[1];
sx q[1];
rz(-0.70901543) q[1];
sx q[1];
rz(0.11911094) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8934694) q[3];
sx q[3];
rz(-2.8827658) q[3];
sx q[3];
rz(2.5138574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36491498) q[2];
sx q[2];
rz(-1.9675576) q[2];
sx q[2];
rz(-0.066369973) q[2];
rz(-2.5868609) q[3];
sx q[3];
rz(-1.4812255) q[3];
sx q[3];
rz(1.5167282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944815) q[0];
sx q[0];
rz(-0.28488657) q[0];
sx q[0];
rz(0.77199212) q[0];
rz(0.66328612) q[1];
sx q[1];
rz(-2.3978077) q[1];
sx q[1];
rz(0.1276806) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049162) q[0];
sx q[0];
rz(-0.66670115) q[0];
sx q[0];
rz(-0.0078860869) q[0];
rz(-pi) q[1];
rz(1.9468609) q[2];
sx q[2];
rz(-2.0916307) q[2];
sx q[2];
rz(-1.2729046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7230715) q[1];
sx q[1];
rz(-1.9772031) q[1];
sx q[1];
rz(2.0416829) q[1];
rz(-0.33416602) q[3];
sx q[3];
rz(-2.6892956) q[3];
sx q[3];
rz(2.3826007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53106236) q[2];
sx q[2];
rz(-0.9610815) q[2];
sx q[2];
rz(-0.5985716) q[2];
rz(-2.5564204) q[3];
sx q[3];
rz(-1.7681311) q[3];
sx q[3];
rz(-0.055805512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7912306) q[0];
sx q[0];
rz(-1.6329916) q[0];
sx q[0];
rz(0.97306657) q[0];
rz(-0.19634253) q[1];
sx q[1];
rz(-1.1390431) q[1];
sx q[1];
rz(1.6729209) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40590433) q[0];
sx q[0];
rz(-2.2668512) q[0];
sx q[0];
rz(0.52175867) q[0];
rz(-pi) q[1];
rz(-2.0765012) q[2];
sx q[2];
rz(-1.7943766) q[2];
sx q[2];
rz(-0.93709842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96983671) q[1];
sx q[1];
rz(-1.7431269) q[1];
sx q[1];
rz(1.8755765) q[1];
rz(-pi) q[2];
rz(-3.0861868) q[3];
sx q[3];
rz(-2.1089777) q[3];
sx q[3];
rz(0.84780592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4520182) q[2];
sx q[2];
rz(-1.3268027) q[2];
sx q[2];
rz(0.16239521) q[2];
rz(-2.0116122) q[3];
sx q[3];
rz(-2.2584848) q[3];
sx q[3];
rz(1.1348772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75055403) q[0];
sx q[0];
rz(-0.51893187) q[0];
sx q[0];
rz(-0.050405141) q[0];
rz(-0.16818908) q[1];
sx q[1];
rz(-1.5805809) q[1];
sx q[1];
rz(-2.2680297) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45437434) q[0];
sx q[0];
rz(-2.4719387) q[0];
sx q[0];
rz(-0.39910103) q[0];
x q[1];
rz(1.4292154) q[2];
sx q[2];
rz(-1.537286) q[2];
sx q[2];
rz(0.17591116) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.37994036) q[1];
sx q[1];
rz(-2.4177175) q[1];
sx q[1];
rz(-0.70751247) q[1];
rz(-pi) q[2];
rz(0.59385517) q[3];
sx q[3];
rz(-2.7866413) q[3];
sx q[3];
rz(-1.6194199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0198589) q[2];
sx q[2];
rz(-0.21616082) q[2];
sx q[2];
rz(-2.7609694) q[2];
rz(-2.5753283) q[3];
sx q[3];
rz(-1.4004613) q[3];
sx q[3];
rz(0.80999058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.64067632) q[0];
sx q[0];
rz(-2.930142) q[0];
sx q[0];
rz(-2.7389615) q[0];
rz(-1.7598049) q[1];
sx q[1];
rz(-2.333162) q[1];
sx q[1];
rz(-0.056338739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005971) q[0];
sx q[0];
rz(-2.6133446) q[0];
sx q[0];
rz(-2.2597367) q[0];
x q[1];
rz(-3.10059) q[2];
sx q[2];
rz(-1.1004538) q[2];
sx q[2];
rz(-0.13861632) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.496324) q[1];
sx q[1];
rz(-0.56425112) q[1];
sx q[1];
rz(-1.5936046) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5272156) q[3];
sx q[3];
rz(-2.3639332) q[3];
sx q[3];
rz(1.7462424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23400433) q[2];
sx q[2];
rz(-1.7981217) q[2];
sx q[2];
rz(-2.6410356) q[2];
rz(0.060700011) q[3];
sx q[3];
rz(-1.3541636) q[3];
sx q[3];
rz(0.48262706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57182264) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(1.7806336) q[0];
rz(0.743615) q[1];
sx q[1];
rz(-1.7230325) q[1];
sx q[1];
rz(-0.13882151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1847938) q[0];
sx q[0];
rz(-1.0652055) q[0];
sx q[0];
rz(0.57539815) q[0];
x q[1];
rz(0.46625579) q[2];
sx q[2];
rz(-1.0276349) q[2];
sx q[2];
rz(-2.0838497) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.946859) q[1];
sx q[1];
rz(-1.7725962) q[1];
sx q[1];
rz(-1.9029593) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0984135) q[3];
sx q[3];
rz(-2.7076027) q[3];
sx q[3];
rz(-1.4066309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4678141) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(-0.71411258) q[2];
rz(0.95947391) q[3];
sx q[3];
rz(-1.4941314) q[3];
sx q[3];
rz(2.1625471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7239083) q[0];
sx q[0];
rz(-0.67643106) q[0];
sx q[0];
rz(-3.0133001) q[0];
rz(-0.3872321) q[1];
sx q[1];
rz(-2.5435244) q[1];
sx q[1];
rz(1.3234352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5538226) q[0];
sx q[0];
rz(-2.9284796) q[0];
sx q[0];
rz(-1.0928667) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5495731) q[2];
sx q[2];
rz(-0.97416234) q[2];
sx q[2];
rz(-0.48812619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8534996) q[1];
sx q[1];
rz(-2.9969831) q[1];
sx q[1];
rz(-1.4734731) q[1];
rz(-pi) q[2];
rz(-1.7701444) q[3];
sx q[3];
rz(-2.2888765) q[3];
sx q[3];
rz(-1.9950641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6748176) q[2];
sx q[2];
rz(-1.6206436) q[2];
sx q[2];
rz(-0.20469323) q[2];
rz(-1.8681059) q[3];
sx q[3];
rz(-2.4114362) q[3];
sx q[3];
rz(1.5654806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765163) q[0];
sx q[0];
rz(-2.5193494) q[0];
sx q[0];
rz(-2.9283071) q[0];
rz(-2.9215096) q[1];
sx q[1];
rz(-1.2525696) q[1];
sx q[1];
rz(1.782104) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3598571) q[0];
sx q[0];
rz(-0.018642519) q[0];
sx q[0];
rz(-2.1117979) q[0];
rz(-1.8259117) q[2];
sx q[2];
rz(-1.1099225) q[2];
sx q[2];
rz(-0.67042353) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9656187) q[1];
sx q[1];
rz(-0.63116628) q[1];
sx q[1];
rz(2.1793773) q[1];
x q[2];
rz(-0.90057365) q[3];
sx q[3];
rz(-1.2069205) q[3];
sx q[3];
rz(-0.28387851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.517211) q[2];
sx q[2];
rz(-1.5456079) q[2];
sx q[2];
rz(0.33779302) q[2];
rz(1.7289915) q[3];
sx q[3];
rz(-1.9905041) q[3];
sx q[3];
rz(0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74078858) q[0];
sx q[0];
rz(-0.82397006) q[0];
sx q[0];
rz(-1.9116221) q[0];
rz(-2.7578655) q[1];
sx q[1];
rz(-1.4839254) q[1];
sx q[1];
rz(0.76212777) q[1];
rz(2.4356213) q[2];
sx q[2];
rz(-2.6298475) q[2];
sx q[2];
rz(2.0012326) q[2];
rz(1.1521793) q[3];
sx q[3];
rz(-2.5910196) q[3];
sx q[3];
rz(-2.2923727) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
