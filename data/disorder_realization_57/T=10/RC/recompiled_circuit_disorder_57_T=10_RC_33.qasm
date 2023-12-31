OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(3.951374) q[0];
sx q[0];
rz(9.9561719) q[0];
rz(-2.9040789) q[1];
sx q[1];
rz(-1.3637435) q[1];
sx q[1];
rz(-1.9385424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42703585) q[0];
sx q[0];
rz(-1.7622951) q[0];
sx q[0];
rz(-1.7869851) q[0];
rz(-pi) q[1];
rz(0.7026303) q[2];
sx q[2];
rz(-0.32877562) q[2];
sx q[2];
rz(-1.6871014) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14129278) q[1];
sx q[1];
rz(-1.80559) q[1];
sx q[1];
rz(1.512212) q[1];
rz(-pi) q[2];
rz(-3.0193127) q[3];
sx q[3];
rz(-0.29502007) q[3];
sx q[3];
rz(2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(-0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99877015) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(-2.3221827) q[0];
rz(-2.8557414) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(1.8751289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0985384) q[0];
sx q[0];
rz(-0.96999723) q[0];
sx q[0];
rz(2.0354802) q[0];
rz(-pi) q[1];
rz(-2.0882294) q[2];
sx q[2];
rz(-2.6839163) q[2];
sx q[2];
rz(-2.6982754) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4283735) q[1];
sx q[1];
rz(-1.1717147) q[1];
sx q[1];
rz(2.748511) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5760848) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(-1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9849898) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.9796237) q[2];
rz(3.0544288) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(-0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6067628) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(-1.3820232) q[1];
sx q[1];
rz(-1.8654114) q[1];
sx q[1];
rz(2.4940925) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.203813) q[0];
sx q[0];
rz(-2.505629) q[0];
sx q[0];
rz(-0.075809191) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0258255) q[2];
sx q[2];
rz(-0.93169824) q[2];
sx q[2];
rz(0.35312032) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6806879) q[1];
sx q[1];
rz(-2.540743) q[1];
sx q[1];
rz(-2.2520573) q[1];
rz(0.87326761) q[3];
sx q[3];
rz(-2.0487361) q[3];
sx q[3];
rz(-2.0877116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(-2.1598699) q[3];
sx q[3];
rz(-1.862062) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(-0.87189829) q[0];
rz(-2.7543228) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(2.8505468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29845333) q[0];
sx q[0];
rz(-1.342134) q[0];
sx q[0];
rz(-2.162621) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8486738) q[2];
sx q[2];
rz(-1.2625164) q[2];
sx q[2];
rz(0.76383797) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17109045) q[1];
sx q[1];
rz(-0.88741747) q[1];
sx q[1];
rz(0.6716397) q[1];
rz(-1.8157186) q[3];
sx q[3];
rz(-2.5513253) q[3];
sx q[3];
rz(-3.0274689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3691833) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(-0.54405653) q[2];
rz(2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0020224) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(-2.496526) q[0];
rz(-0.6257261) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(2.5114139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7932927) q[0];
sx q[0];
rz(-0.59069809) q[0];
sx q[0];
rz(2.9212055) q[0];
rz(-pi) q[1];
rz(-0.39406392) q[2];
sx q[2];
rz(-2.2820018) q[2];
sx q[2];
rz(-1.994339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1944151) q[1];
sx q[1];
rz(-1.0395323) q[1];
sx q[1];
rz(-2.8814425) q[1];
rz(-pi) q[2];
rz(2.9003521) q[3];
sx q[3];
rz(-0.57096982) q[3];
sx q[3];
rz(2.8591448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3877635) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(-1.7774263) q[2];
rz(-1.2498614) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867778) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-0.71682799) q[0];
rz(-2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(0.81370083) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98052927) q[0];
sx q[0];
rz(-1.6335532) q[0];
sx q[0];
rz(-1.784523) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5021624) q[2];
sx q[2];
rz(-0.8317906) q[2];
sx q[2];
rz(2.8384428) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8061714) q[1];
sx q[1];
rz(-0.85039447) q[1];
sx q[1];
rz(2.6909188) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86088647) q[3];
sx q[3];
rz(-2.6885899) q[3];
sx q[3];
rz(2.3080254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8852691) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(1.1876594) q[2];
rz(0.93196431) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-0.37117547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9880144) q[0];
sx q[0];
rz(-2.1126641) q[0];
sx q[0];
rz(-2.4555092) q[0];
rz(-1.6185121) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(-0.66326052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30090573) q[0];
sx q[0];
rz(-1.5468883) q[0];
sx q[0];
rz(2.2211214) q[0];
rz(-pi) q[1];
rz(2.8510423) q[2];
sx q[2];
rz(-2.4834589) q[2];
sx q[2];
rz(-1.5302637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7546332) q[1];
sx q[1];
rz(-2.6135751) q[1];
sx q[1];
rz(-1.1623043) q[1];
rz(-pi) q[2];
rz(0.56179629) q[3];
sx q[3];
rz(-1.1601163) q[3];
sx q[3];
rz(0.79515776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.474581) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(2.1298501) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(-2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(-3.0342039) q[0];
rz(-0.30474162) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-2.6838141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4360355) q[0];
sx q[0];
rz(-1.2495263) q[0];
sx q[0];
rz(0.9106439) q[0];
rz(0.7465676) q[2];
sx q[2];
rz(-1.1636359) q[2];
sx q[2];
rz(-2.6395869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.01955186) q[1];
sx q[1];
rz(-1.3380088) q[1];
sx q[1];
rz(-1.9368534) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55191603) q[3];
sx q[3];
rz(-1.8822) q[3];
sx q[3];
rz(0.50406721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36499873) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.7101074) q[2];
rz(1.4252023) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(-1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872221) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(-0.95348683) q[0];
rz(1.3257239) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(0.58473933) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23823243) q[0];
sx q[0];
rz(-1.0833217) q[0];
sx q[0];
rz(1.4438629) q[0];
x q[1];
rz(0.54359162) q[2];
sx q[2];
rz(-1.2541176) q[2];
sx q[2];
rz(-1.8900332) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7331839) q[1];
sx q[1];
rz(-1.5015263) q[1];
sx q[1];
rz(-0.44002156) q[1];
x q[2];
rz(-1.8381848) q[3];
sx q[3];
rz(-1.7087473) q[3];
sx q[3];
rz(1.4064058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(-2.8653223) q[2];
rz(-0.55109465) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(0.38366693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0722512) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(-1.9845225) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(-1.5225333) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9690543) q[0];
sx q[0];
rz(-1.8803055) q[0];
sx q[0];
rz(-2.4597416) q[0];
x q[1];
rz(1.2298146) q[2];
sx q[2];
rz(-2.2298862) q[2];
sx q[2];
rz(-1.9552719) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9177502) q[1];
sx q[1];
rz(-0.94864861) q[1];
sx q[1];
rz(1.7734581) q[1];
rz(-pi) q[2];
rz(-0.75928648) q[3];
sx q[3];
rz(-2.5873103) q[3];
sx q[3];
rz(-2.5169773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(-0.84214169) q[2];
rz(1.5367674) q[3];
sx q[3];
rz(-1.3804881) q[3];
sx q[3];
rz(3.1150637) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90606541) q[0];
sx q[0];
rz(-1.9530095) q[0];
sx q[0];
rz(-0.45146913) q[0];
rz(0.72262598) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(1.4459544) q[2];
sx q[2];
rz(-0.21077934) q[2];
sx q[2];
rz(-1.6817033) q[2];
rz(0.87993965) q[3];
sx q[3];
rz(-1.5094681) q[3];
sx q[3];
rz(-1.5916011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
