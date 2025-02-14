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
rz(-1.3396076) q[0];
sx q[0];
rz(-0.34914246) q[0];
sx q[0];
rz(2.1139297) q[0];
rz(-0.034962058) q[1];
sx q[1];
rz(-1.0171913) q[1];
sx q[1];
rz(-1.4453759) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6562614) q[0];
sx q[0];
rz(-2.4527271) q[0];
sx q[0];
rz(0.553225) q[0];
x q[1];
rz(0.18741636) q[2];
sx q[2];
rz(-0.18067154) q[2];
sx q[2];
rz(-0.60127163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0398324) q[1];
sx q[1];
rz(-1.2451485) q[1];
sx q[1];
rz(-1.7937307) q[1];
x q[2];
rz(-2.7463116) q[3];
sx q[3];
rz(-2.3594405) q[3];
sx q[3];
rz(-2.8349193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8029636) q[2];
sx q[2];
rz(-2.6821319) q[2];
sx q[2];
rz(-2.6098693) q[2];
rz(-1.7470597) q[3];
sx q[3];
rz(-1.2732384) q[3];
sx q[3];
rz(-1.9664221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8748473) q[0];
sx q[0];
rz(-1.5044455) q[0];
sx q[0];
rz(-0.7919842) q[0];
rz(-2.2199471) q[1];
sx q[1];
rz(-1.6519203) q[1];
sx q[1];
rz(1.1080866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32852325) q[0];
sx q[0];
rz(-1.4591154) q[0];
sx q[0];
rz(-2.172676) q[0];
rz(-pi) q[1];
rz(-0.048594995) q[2];
sx q[2];
rz(-0.94779769) q[2];
sx q[2];
rz(0.6374661) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11985677) q[1];
sx q[1];
rz(-0.83470063) q[1];
sx q[1];
rz(-1.0480919) q[1];
x q[2];
rz(-1.7407623) q[3];
sx q[3];
rz(-2.8208591) q[3];
sx q[3];
rz(-0.50839822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7532928) q[2];
sx q[2];
rz(-1.029195) q[2];
sx q[2];
rz(-1.8886867) q[2];
rz(-0.62475723) q[3];
sx q[3];
rz(-1.8936936) q[3];
sx q[3];
rz(-2.6784082) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4598292) q[0];
sx q[0];
rz(-0.1777996) q[0];
sx q[0];
rz(2.8681927) q[0];
rz(3.0699442) q[1];
sx q[1];
rz(-2.007808) q[1];
sx q[1];
rz(1.6260446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1119627) q[0];
sx q[0];
rz(-0.66034895) q[0];
sx q[0];
rz(-0.65159722) q[0];
rz(-2.5070174) q[2];
sx q[2];
rz(-0.56152841) q[2];
sx q[2];
rz(0.47088366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7715986) q[1];
sx q[1];
rz(-0.99111667) q[1];
sx q[1];
rz(1.3808151) q[1];
rz(-pi) q[2];
rz(-1.9380366) q[3];
sx q[3];
rz(-0.62588309) q[3];
sx q[3];
rz(0.82928951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7800954) q[2];
sx q[2];
rz(-0.7146892) q[2];
sx q[2];
rz(-1.7717465) q[2];
rz(2.0731549) q[3];
sx q[3];
rz(-1.3911824) q[3];
sx q[3];
rz(0.95503241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84900981) q[0];
sx q[0];
rz(-2.7136901) q[0];
sx q[0];
rz(0.12246116) q[0];
rz(0.017008688) q[1];
sx q[1];
rz(-1.5127134) q[1];
sx q[1];
rz(0.88517991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1830782) q[0];
sx q[0];
rz(-1.7066656) q[0];
sx q[0];
rz(2.0625173) q[0];
x q[1];
rz(-0.34220747) q[2];
sx q[2];
rz(-0.86538431) q[2];
sx q[2];
rz(-0.17555412) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7834691) q[1];
sx q[1];
rz(-1.1162288) q[1];
sx q[1];
rz(-0.44444167) q[1];
rz(-pi) q[2];
rz(-1.3766798) q[3];
sx q[3];
rz(-2.2633865) q[3];
sx q[3];
rz(3.0339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.433832) q[2];
sx q[2];
rz(-1.8261352) q[2];
sx q[2];
rz(0.07746499) q[2];
rz(-0.42099434) q[3];
sx q[3];
rz(-1.9546031) q[3];
sx q[3];
rz(0.25119701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.0635327) q[0];
sx q[0];
rz(-3.1227626) q[0];
sx q[0];
rz(2.4596762) q[0];
rz(-0.49184999) q[1];
sx q[1];
rz(-1.0268772) q[1];
sx q[1];
rz(1.1600201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60940015) q[0];
sx q[0];
rz(-0.74157292) q[0];
sx q[0];
rz(-1.8511008) q[0];
rz(-0.087988532) q[2];
sx q[2];
rz(-0.57214979) q[2];
sx q[2];
rz(2.3038626) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.034879) q[1];
sx q[1];
rz(-2.1886161) q[1];
sx q[1];
rz(0.21150963) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8468708) q[3];
sx q[3];
rz(-1.975276) q[3];
sx q[3];
rz(-1.7691139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68449768) q[2];
sx q[2];
rz(-2.2726077) q[2];
sx q[2];
rz(-1.0082461) q[2];
rz(-1.6393939) q[3];
sx q[3];
rz(-1.1049756) q[3];
sx q[3];
rz(2.8777299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0428001) q[0];
sx q[0];
rz(-1.5941987) q[0];
sx q[0];
rz(0.40476009) q[0];
rz(-2.3233991) q[1];
sx q[1];
rz(-1.300756) q[1];
sx q[1];
rz(0.71664804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.19159) q[0];
sx q[0];
rz(-1.5609589) q[0];
sx q[0];
rz(0.27920009) q[0];
x q[1];
rz(0.54714142) q[2];
sx q[2];
rz(-0.84274492) q[2];
sx q[2];
rz(2.0942618) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7886519) q[1];
sx q[1];
rz(-1.8308906) q[1];
sx q[1];
rz(-2.6648894) q[1];
rz(-pi) q[2];
rz(0.53933177) q[3];
sx q[3];
rz(-1.2507273) q[3];
sx q[3];
rz(-2.0254997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8019668) q[2];
sx q[2];
rz(-2.5012987) q[2];
sx q[2];
rz(0.65529811) q[2];
rz(-2.7845434) q[3];
sx q[3];
rz(-1.7506295) q[3];
sx q[3];
rz(2.6220139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.201467) q[0];
sx q[0];
rz(-2.8253912) q[0];
sx q[0];
rz(-2.1215718) q[0];
rz(0.54234281) q[1];
sx q[1];
rz(-1.1154563) q[1];
sx q[1];
rz(-1.7592336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64847222) q[0];
sx q[0];
rz(-2.1628404) q[0];
sx q[0];
rz(-0.73541321) q[0];
rz(-pi) q[1];
rz(0.57047259) q[2];
sx q[2];
rz(-2.8920724) q[2];
sx q[2];
rz(-0.0970627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51782896) q[1];
sx q[1];
rz(-1.4925071) q[1];
sx q[1];
rz(-1.109156) q[1];
rz(-pi) q[2];
rz(1.5685308) q[3];
sx q[3];
rz(-1.6544154) q[3];
sx q[3];
rz(2.9155242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0688087) q[2];
sx q[2];
rz(-1.5280318) q[2];
sx q[2];
rz(2.7911348) q[2];
rz(2.6751878) q[3];
sx q[3];
rz(-0.91752183) q[3];
sx q[3];
rz(-0.94356999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7826409) q[0];
sx q[0];
rz(-2.751001) q[0];
sx q[0];
rz(-2.1851831) q[0];
rz(1.4495173) q[1];
sx q[1];
rz(-2.3608975) q[1];
sx q[1];
rz(-0.67109674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.886837) q[0];
sx q[0];
rz(-0.16567812) q[0];
sx q[0];
rz(-2.1149693) q[0];
rz(-pi) q[1];
rz(-0.32197774) q[2];
sx q[2];
rz(-2.5484332) q[2];
sx q[2];
rz(-0.71074394) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33167142) q[1];
sx q[1];
rz(-1.6345895) q[1];
sx q[1];
rz(-1.2170736) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28658925) q[3];
sx q[3];
rz(-2.4169888) q[3];
sx q[3];
rz(-2.2514908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2676919) q[2];
sx q[2];
rz(-2.716422) q[2];
sx q[2];
rz(-2.0153866) q[2];
rz(-2.8113484) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(-0.11708524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983593) q[0];
sx q[0];
rz(-0.57906228) q[0];
sx q[0];
rz(0.057057127) q[0];
rz(-0.13380274) q[1];
sx q[1];
rz(-2.0963142) q[1];
sx q[1];
rz(-0.71279508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4445515) q[0];
sx q[0];
rz(-1.6336055) q[0];
sx q[0];
rz(2.8256326) q[0];
rz(-pi) q[1];
rz(-0.78929094) q[2];
sx q[2];
rz(-0.54614353) q[2];
sx q[2];
rz(-0.30420732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1755029) q[1];
sx q[1];
rz(-1.1017297) q[1];
sx q[1];
rz(1.1203946) q[1];
rz(-pi) q[2];
rz(-2.2650293) q[3];
sx q[3];
rz(-1.6329771) q[3];
sx q[3];
rz(-0.93557318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.628525) q[2];
sx q[2];
rz(-1.2863938) q[2];
sx q[2];
rz(0.081550278) q[2];
rz(1.2201803) q[3];
sx q[3];
rz(-0.27118513) q[3];
sx q[3];
rz(2.0512106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9995025) q[0];
sx q[0];
rz(-1.6195848) q[0];
sx q[0];
rz(-1.581544) q[0];
rz(2.0060495) q[1];
sx q[1];
rz(-1.138843) q[1];
sx q[1];
rz(1.1467689) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97871507) q[0];
sx q[0];
rz(-2.204127) q[0];
sx q[0];
rz(-0.94453728) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76374526) q[2];
sx q[2];
rz(-0.95852214) q[2];
sx q[2];
rz(2.6291922) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9830977) q[1];
sx q[1];
rz(-1.4870411) q[1];
sx q[1];
rz(2.9309446) q[1];
x q[2];
rz(1.0057529) q[3];
sx q[3];
rz(-1.4970137) q[3];
sx q[3];
rz(-2.0948054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59209383) q[2];
sx q[2];
rz(-1.3201951) q[2];
sx q[2];
rz(2.7756694) q[2];
rz(-2.7538815) q[3];
sx q[3];
rz(-1.1185027) q[3];
sx q[3];
rz(-1.6754735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53030071) q[0];
sx q[0];
rz(-2.3983751) q[0];
sx q[0];
rz(-0.30839738) q[0];
rz(-2.3920234) q[1];
sx q[1];
rz(-2.0105965) q[1];
sx q[1];
rz(-1.2731193) q[1];
rz(0.34229924) q[2];
sx q[2];
rz(-0.93955775) q[2];
sx q[2];
rz(-1.0472236) q[2];
rz(0.6324296) q[3];
sx q[3];
rz(-1.8550183) q[3];
sx q[3];
rz(0.99799533) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
