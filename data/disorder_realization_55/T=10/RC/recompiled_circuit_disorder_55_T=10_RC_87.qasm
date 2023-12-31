OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(4.6580553) q[0];
sx q[0];
rz(9.1604995) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(5.073054) q[1];
sx q[1];
rz(10.160025) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7289294) q[0];
sx q[0];
rz(-0.67617765) q[0];
sx q[0];
rz(2.9039608) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5392883) q[2];
sx q[2];
rz(-0.77568433) q[2];
sx q[2];
rz(1.3210981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9303317) q[1];
sx q[1];
rz(-0.32713612) q[1];
sx q[1];
rz(1.0717908) q[1];
x q[2];
rz(1.5427038) q[3];
sx q[3];
rz(-0.99053226) q[3];
sx q[3];
rz(2.7147646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66951093) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(-0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(-2.7052178) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-2.8754821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9929745) q[0];
sx q[0];
rz(-2.2668112) q[0];
sx q[0];
rz(-2.6398185) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9782148) q[2];
sx q[2];
rz(-1.0014357) q[2];
sx q[2];
rz(-0.53158224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38763603) q[1];
sx q[1];
rz(-1.0322744) q[1];
sx q[1];
rz(2.3073879) q[1];
x q[2];
rz(0.9640785) q[3];
sx q[3];
rz(-1.7214516) q[3];
sx q[3];
rz(1.7863303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(0.51149386) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(-2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(0.92873746) q[0];
rz(1.7354895) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.3471289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0179694) q[0];
sx q[0];
rz(-2.6051913) q[0];
sx q[0];
rz(-0.53039741) q[0];
rz(1.2855661) q[2];
sx q[2];
rz(-1.4189548) q[2];
sx q[2];
rz(1.0361995) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.320142) q[1];
sx q[1];
rz(-2.4648033) q[1];
sx q[1];
rz(2.0435145) q[1];
x q[2];
rz(-2.9903528) q[3];
sx q[3];
rz(-0.86942196) q[3];
sx q[3];
rz(2.410694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1469664) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(-1.6195126) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(-0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(2.175892) q[0];
rz(0.72215885) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(0.55975634) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7402732) q[0];
sx q[0];
rz(-1.4896605) q[0];
sx q[0];
rz(-0.42663891) q[0];
rz(-0.83300029) q[2];
sx q[2];
rz(-2.2350395) q[2];
sx q[2];
rz(-0.26495648) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30448118) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(-2.1865305) q[1];
x q[2];
rz(1.8074606) q[3];
sx q[3];
rz(-2.3458614) q[3];
sx q[3];
rz(-2.3468897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7136148) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(-2.8811841) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.150862) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(2.7752303) q[0];
rz(-1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(-2.7979134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8415547) q[0];
sx q[0];
rz(-0.90111387) q[0];
sx q[0];
rz(0.71398736) q[0];
rz(-0.57858606) q[2];
sx q[2];
rz(-0.95539504) q[2];
sx q[2];
rz(1.3784642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3949563) q[1];
sx q[1];
rz(-2.6379105) q[1];
sx q[1];
rz(1.7240745) q[1];
rz(-pi) q[2];
x q[2];
rz(2.214659) q[3];
sx q[3];
rz(-1.5537795) q[3];
sx q[3];
rz(3.0683558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7498103) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(-3.0878477) q[2];
rz(1.7371477) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(0.25973928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8206257) q[0];
sx q[0];
rz(-0.48829406) q[0];
sx q[0];
rz(-2.2886306) q[0];
rz(-pi) q[1];
rz(1.514228) q[2];
sx q[2];
rz(-2.8048189) q[2];
sx q[2];
rz(2.3677504) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7139587) q[1];
sx q[1];
rz(-2.7379588) q[1];
sx q[1];
rz(-0.9637109) q[1];
rz(-pi) q[2];
rz(1.5335347) q[3];
sx q[3];
rz(-0.44964368) q[3];
sx q[3];
rz(-0.26526181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(-3.0409813) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(-1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1473734) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(-2.8314262) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(-2.535634) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6771686) q[0];
sx q[0];
rz(-2.4806528) q[0];
sx q[0];
rz(-2.7695157) q[0];
x q[1];
rz(1.7714959) q[2];
sx q[2];
rz(-1.8569274) q[2];
sx q[2];
rz(-1.2933033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0026605) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(2.1041811) q[1];
rz(-pi) q[2];
rz(2.1171655) q[3];
sx q[3];
rz(-1.1064648) q[3];
sx q[3];
rz(1.6950316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24017748) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(-2.288738) q[2];
rz(-1.3700221) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119499) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(-1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(-3.0775552) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56211573) q[0];
sx q[0];
rz(-2.1770283) q[0];
sx q[0];
rz(2.9129145) q[0];
rz(0.30631752) q[2];
sx q[2];
rz(-2.0567354) q[2];
sx q[2];
rz(1.475032) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36754164) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(2.7864085) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1021032) q[3];
sx q[3];
rz(-1.8825304) q[3];
sx q[3];
rz(1.5962275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.091207592) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-2.2122038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9973307) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147946) q[0];
sx q[0];
rz(-0.1495805) q[0];
sx q[0];
rz(3.0339255) q[0];
x q[1];
rz(2.1986507) q[2];
sx q[2];
rz(-1.1254416) q[2];
sx q[2];
rz(3.359059e-05) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7334426) q[1];
sx q[1];
rz(-2.8996455) q[1];
sx q[1];
rz(1.5223632) q[1];
rz(2.7752152) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(0.60929326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.613712) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(0.17364994) q[2];
rz(2.8052143) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(-2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.4655112) q[0];
rz(-2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-2.5691659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51047072) q[0];
sx q[0];
rz(-2.0513751) q[0];
sx q[0];
rz(-2.5170588) q[0];
rz(-1.1119214) q[2];
sx q[2];
rz(-1.8177114) q[2];
sx q[2];
rz(-1.0690451) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.010667) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(0.24433498) q[1];
x q[2];
rz(-2.8045373) q[3];
sx q[3];
rz(-2.4042077) q[3];
sx q[3];
rz(2.9205703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77999014) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(0.53722107) q[2];
rz(-2.0843263) q[3];
sx q[3];
rz(-2.2500762) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.854241) q[0];
sx q[0];
rz(-1.9762522) q[0];
sx q[0];
rz(1.5594788) q[0];
rz(-2.6782425) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(-2.2255185) q[2];
sx q[2];
rz(-1.6486042) q[2];
sx q[2];
rz(-0.83124607) q[2];
rz(-1.6871917) q[3];
sx q[3];
rz(-2.6506861) q[3];
sx q[3];
rz(-0.41674137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
