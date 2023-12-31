OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6131634) q[0];
sx q[0];
rz(-2.0818721) q[0];
sx q[0];
rz(-0.73097316) q[0];
rz(1.641474) q[1];
sx q[1];
rz(-1.0348231) q[1];
sx q[1];
rz(2.1980481) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7109414) q[0];
sx q[0];
rz(-1.5358155) q[0];
sx q[0];
rz(3.0745688) q[0];
rz(-pi) q[1];
rz(-3.0303454) q[2];
sx q[2];
rz(-2.4876378) q[2];
sx q[2];
rz(-2.7788869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2343443) q[1];
sx q[1];
rz(-0.53011361) q[1];
sx q[1];
rz(-2.2905473) q[1];
rz(2.5991873) q[3];
sx q[3];
rz(-0.35202682) q[3];
sx q[3];
rz(2.5792518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(-1.250766) q[2];
rz(1.7154153) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(-2.1616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0086867) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(2.5426478) q[0];
rz(1.3409746) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(0.96639955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7705298) q[0];
sx q[0];
rz(-2.3822228) q[0];
sx q[0];
rz(0.66803996) q[0];
x q[1];
rz(-2.5194089) q[2];
sx q[2];
rz(-1.5303648) q[2];
sx q[2];
rz(1.7379023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2981373) q[1];
sx q[1];
rz(-1.6846488) q[1];
sx q[1];
rz(-1.650412) q[1];
x q[2];
rz(-0.45774777) q[3];
sx q[3];
rz(-0.80349892) q[3];
sx q[3];
rz(-0.97704923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(-1.1931233) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10467228) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(1.954129) q[0];
rz(1.9056412) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.3175861) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0061958) q[0];
sx q[0];
rz(-0.72512308) q[0];
sx q[0];
rz(2.8185185) q[0];
rz(-pi) q[1];
x q[1];
rz(2.927711) q[2];
sx q[2];
rz(-1.7638532) q[2];
sx q[2];
rz(-2.4436827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31955645) q[1];
sx q[1];
rz(-1.6720547) q[1];
sx q[1];
rz(-0.42145573) q[1];
x q[2];
rz(1.0533603) q[3];
sx q[3];
rz(-1.6728757) q[3];
sx q[3];
rz(1.3641588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(0.2066361) q[2];
rz(0.7080428) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77465039) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(-0.88622093) q[0];
rz(1.0097424) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(1.2264235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013989) q[0];
sx q[0];
rz(-1.6370019) q[0];
sx q[0];
rz(2.2982236) q[0];
rz(-pi) q[1];
rz(2.6817276) q[2];
sx q[2];
rz(-2.2042639) q[2];
sx q[2];
rz(1.2666653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1434506) q[1];
sx q[1];
rz(-1.3610024) q[1];
sx q[1];
rz(-1.3694805) q[1];
x q[2];
rz(0.48638101) q[3];
sx q[3];
rz(-1.2380621) q[3];
sx q[3];
rz(0.90852028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(-0.99299661) q[2];
rz(1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(2.5591154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5486149) q[0];
sx q[0];
rz(-0.60615221) q[0];
sx q[0];
rz(-1.4447601) q[0];
x q[1];
rz(-2.2485579) q[2];
sx q[2];
rz(-1.4544011) q[2];
sx q[2];
rz(-1.2410156) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1867265) q[1];
sx q[1];
rz(-2.426882) q[1];
sx q[1];
rz(3.1123118) q[1];
rz(-pi) q[2];
rz(-1.5037687) q[3];
sx q[3];
rz(-1.9762632) q[3];
sx q[3];
rz(1.0073347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65486583) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(-0.47406667) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24494568) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(-3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(2.0369464) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4545126) q[0];
sx q[0];
rz(-2.2283163) q[0];
sx q[0];
rz(-1.3643273) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32128895) q[2];
sx q[2];
rz(-0.66547223) q[2];
sx q[2];
rz(1.9208391) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.066597477) q[1];
sx q[1];
rz(-2.1146333) q[1];
sx q[1];
rz(2.4165513) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8940582) q[3];
sx q[3];
rz(-2.6544016) q[3];
sx q[3];
rz(0.67684735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8905939) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(0.55523038) q[2];
rz(2.9688719) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5884488) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-2.0297091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33169532) q[0];
sx q[0];
rz(-0.90242093) q[0];
sx q[0];
rz(-2.6033127) q[0];
x q[1];
rz(1.7467473) q[2];
sx q[2];
rz(-3.0627652) q[2];
sx q[2];
rz(-1.4836756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8916787) q[1];
sx q[1];
rz(-2.2644682) q[1];
sx q[1];
rz(1.9338884) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8133409) q[3];
sx q[3];
rz(-1.8040856) q[3];
sx q[3];
rz(-2.1094028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(0.81364441) q[2];
rz(1.404473) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62548816) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.7161436) q[0];
rz(-1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(2.5040748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2256048) q[0];
sx q[0];
rz(-2.2921352) q[0];
sx q[0];
rz(2.3774873) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9279187) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(1.9879607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9962822) q[1];
sx q[1];
rz(-1.1728371) q[1];
sx q[1];
rz(-2.5064962) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5708837) q[3];
sx q[3];
rz(-1.4548886) q[3];
sx q[3];
rz(-1.4012208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0104388) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(-2.0054224) q[2];
rz(1.6561967) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(-1.2095399) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5159601) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(2.4654454) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069106) q[0];
sx q[0];
rz(-2.7043531) q[0];
sx q[0];
rz(-2.8291563) q[0];
rz(-pi) q[1];
rz(0.65281547) q[2];
sx q[2];
rz(-1.6932994) q[2];
sx q[2];
rz(-1.1268238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91190126) q[1];
sx q[1];
rz(-1.9209314) q[1];
sx q[1];
rz(-1.9446816) q[1];
rz(-pi) q[2];
rz(-1.2328524) q[3];
sx q[3];
rz(-0.83331185) q[3];
sx q[3];
rz(2.1144707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(0.96735111) q[2];
rz(1.597065) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230781) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.7369695) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.343095) q[0];
sx q[0];
rz(-1.4416749) q[0];
sx q[0];
rz(-1.9341747) q[0];
rz(-pi) q[1];
x q[1];
rz(2.769906) q[2];
sx q[2];
rz(-2.8166397) q[2];
sx q[2];
rz(1.6290968) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9060668) q[1];
sx q[1];
rz(-1.1065673) q[1];
sx q[1];
rz(-0.81495754) q[1];
x q[2];
rz(2.9910562) q[3];
sx q[3];
rz(-1.5794465) q[3];
sx q[3];
rz(-1.0516402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(2.7434769) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.1595935) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(-1.8691312) q[2];
sx q[2];
rz(-1.3186426) q[2];
sx q[2];
rz(-1.1124055) q[2];
rz(-2.0648099) q[3];
sx q[3];
rz(-1.8643338) q[3];
sx q[3];
rz(-2.1716933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
