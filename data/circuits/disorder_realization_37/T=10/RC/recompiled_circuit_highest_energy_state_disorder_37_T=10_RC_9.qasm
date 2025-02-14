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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(0.45912826) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(4.876457) q[1];
sx q[1];
rz(9.193037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.988294) q[0];
sx q[0];
rz(-0.35044119) q[0];
sx q[0];
rz(-2.9399583) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54643537) q[2];
sx q[2];
rz(-1.6197512) q[2];
sx q[2];
rz(3.0374683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5197088) q[1];
sx q[1];
rz(-2.1823898) q[1];
sx q[1];
rz(2.0382463) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0240293) q[3];
sx q[3];
rz(-0.34053482) q[3];
sx q[3];
rz(-1.3887608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0509402) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(-1.2789307) q[2];
rz(2.2517962) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58562529) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(-2.2157748) q[0];
rz(-2.0766808) q[1];
sx q[1];
rz(-1.5771259) q[1];
sx q[1];
rz(-0.085478641) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8259698) q[0];
sx q[0];
rz(-1.1809826) q[0];
sx q[0];
rz(1.5672383) q[0];
rz(-pi) q[1];
rz(1.0930834) q[2];
sx q[2];
rz(-2.4730549) q[2];
sx q[2];
rz(-2.3195409) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.46175112) q[1];
sx q[1];
rz(-2.0978758) q[1];
sx q[1];
rz(-2.8746469) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0382337) q[3];
sx q[3];
rz(-2.2027443) q[3];
sx q[3];
rz(-1.2480717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3257137) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(-0.386664) q[2];
rz(-1.2169085) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(-2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3306408) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(0.49705848) q[0];
rz(-2.0992384) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(2.846948) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93268585) q[0];
sx q[0];
rz(-2.8021376) q[0];
sx q[0];
rz(-0.20821388) q[0];
rz(-pi) q[1];
rz(1.7435837) q[2];
sx q[2];
rz(-2.2026988) q[2];
sx q[2];
rz(1.8586707) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68874796) q[1];
sx q[1];
rz(-1.1383245) q[1];
sx q[1];
rz(0.8853064) q[1];
x q[2];
rz(1.3623275) q[3];
sx q[3];
rz(-2.5304171) q[3];
sx q[3];
rz(1.1217505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3543388) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(-1.5617237) q[2];
rz(-2.0653557) q[3];
sx q[3];
rz(-1.302224) q[3];
sx q[3];
rz(2.2126183) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314986) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(0.22931799) q[0];
rz(1.6353105) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(-3.07952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6825166) q[0];
sx q[0];
rz(-2.354265) q[0];
sx q[0];
rz(2.3045901) q[0];
rz(-pi) q[1];
rz(-1.1787291) q[2];
sx q[2];
rz(-1.3118366) q[2];
sx q[2];
rz(-2.1990537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0921539) q[1];
sx q[1];
rz(-2.2518932) q[1];
sx q[1];
rz(-1.9309631) q[1];
rz(-pi) q[2];
rz(-1.6639978) q[3];
sx q[3];
rz(-1.3706932) q[3];
sx q[3];
rz(-2.4148108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.092992358) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(-1.3108866) q[2];
rz(0.038289573) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(2.766975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6467317) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(-0.89299655) q[0];
rz(2.112174) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(-3.0457048) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9150118) q[0];
sx q[0];
rz(-1.1553197) q[0];
sx q[0];
rz(-1.0092495) q[0];
rz(-0.070438373) q[2];
sx q[2];
rz(-1.4820921) q[2];
sx q[2];
rz(-0.34102893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3600625) q[1];
sx q[1];
rz(-1.2075829) q[1];
sx q[1];
rz(-2.5468154) q[1];
x q[2];
rz(0.51561098) q[3];
sx q[3];
rz(-1.8985975) q[3];
sx q[3];
rz(2.9189381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15478495) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(-0.42243877) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(-2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7285889) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(3.0244306) q[0];
rz(-1.6475742) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(2.07043) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2334918) q[0];
sx q[0];
rz(-1.0885493) q[0];
sx q[0];
rz(1.890475) q[0];
rz(-0.54212924) q[2];
sx q[2];
rz(-1.8111992) q[2];
sx q[2];
rz(-1.8612651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2177762) q[1];
sx q[1];
rz(-0.72715302) q[1];
sx q[1];
rz(-1.9132861) q[1];
rz(-1.4070687) q[3];
sx q[3];
rz(-1.8843972) q[3];
sx q[3];
rz(-1.1160451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5693207) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(3.102109) q[2];
rz(-3.0651921) q[3];
sx q[3];
rz(-1.971784) q[3];
sx q[3];
rz(-2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70306784) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(1.1514661) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(-1.3075525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9614544) q[0];
sx q[0];
rz(-0.39648065) q[0];
sx q[0];
rz(2.1679818) q[0];
rz(-pi) q[1];
rz(1.6436395) q[2];
sx q[2];
rz(-1.9460287) q[2];
sx q[2];
rz(-0.90922395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89245172) q[1];
sx q[1];
rz(-2.0607407) q[1];
sx q[1];
rz(-0.23995121) q[1];
x q[2];
rz(-2.1767307) q[3];
sx q[3];
rz(-1.6814702) q[3];
sx q[3];
rz(-1.0718653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3211956) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(2.5879228) q[2];
rz(-2.4456444) q[3];
sx q[3];
rz(-1.4724933) q[3];
sx q[3];
rz(1.6863916) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11947908) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(-2.8346862) q[0];
rz(-1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(1.2409522) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4590603) q[0];
sx q[0];
rz(-2.0143565) q[0];
sx q[0];
rz(0.17291594) q[0];
rz(1.9996793) q[2];
sx q[2];
rz(-2.1253573) q[2];
sx q[2];
rz(2.6311324) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7684264) q[1];
sx q[1];
rz(-1.2618999) q[1];
sx q[1];
rz(-1.7684446) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1697024) q[3];
sx q[3];
rz(-2.1425284) q[3];
sx q[3];
rz(-0.16756646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79157311) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(1.8512858) q[2];
rz(2.4604515) q[3];
sx q[3];
rz(-0.9001503) q[3];
sx q[3];
rz(-1.5017989) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8660368) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(1.9679029) q[0];
rz(-0.58700079) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(-3.0618844) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22686887) q[0];
sx q[0];
rz(-0.43926111) q[0];
sx q[0];
rz(1.0160116) q[0];
rz(-pi) q[1];
rz(-1.4585895) q[2];
sx q[2];
rz(-1.7908899) q[2];
sx q[2];
rz(0.28973636) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.798493) q[1];
sx q[1];
rz(-2.1329555) q[1];
sx q[1];
rz(2.4141099) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2740785) q[3];
sx q[3];
rz(-0.78700262) q[3];
sx q[3];
rz(0.84144652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54032636) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(1.3339174) q[2];
rz(1.5909083) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(-2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48225668) q[0];
sx q[0];
rz(-1.5784669) q[0];
sx q[0];
rz(1.2905066) q[0];
rz(0.036529649) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(-1.0443002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.999015) q[0];
sx q[0];
rz(-1.2730838) q[0];
sx q[0];
rz(2.0172869) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9007334) q[2];
sx q[2];
rz(-1.3173027) q[2];
sx q[2];
rz(-1.1546749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.83936939) q[1];
sx q[1];
rz(-1.4166792) q[1];
sx q[1];
rz(-0.83854143) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.812522) q[3];
sx q[3];
rz(-1.6794723) q[3];
sx q[3];
rz(-1.4377126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2269939) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(1.8168137) q[2];
rz(-2.3850208) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751547) q[0];
sx q[0];
rz(-1.403724) q[0];
sx q[0];
rz(-2.4028461) q[0];
rz(1.4229763) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(1.1463317) q[2];
sx q[2];
rz(-2.2590315) q[2];
sx q[2];
rz(-3.0207241) q[2];
rz(2.041009) q[3];
sx q[3];
rz(-2.1775424) q[3];
sx q[3];
rz(1.2829124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
