OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72803175) q[0];
sx q[0];
rz(-0.95946884) q[0];
sx q[0];
rz(2.5897107) q[0];
rz(2.7565487) q[1];
sx q[1];
rz(-1.7793964) q[1];
sx q[1];
rz(2.722932) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4872348) q[0];
sx q[0];
rz(-2.4208768) q[0];
sx q[0];
rz(1.1281752) q[0];
x q[1];
rz(-0.35595591) q[2];
sx q[2];
rz(-2.666143) q[2];
sx q[2];
rz(-0.33671492) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.31778279) q[1];
sx q[1];
rz(-1.6595073) q[1];
sx q[1];
rz(0.20388468) q[1];
x q[2];
rz(-1.5403662) q[3];
sx q[3];
rz(-1.55369) q[3];
sx q[3];
rz(-0.78029437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.60575134) q[2];
sx q[2];
rz(-1.302482) q[2];
sx q[2];
rz(-0.38899404) q[2];
rz(-2.244106) q[3];
sx q[3];
rz(-2.5806081) q[3];
sx q[3];
rz(1.8628023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9101343) q[0];
sx q[0];
rz(-2.1899905) q[0];
sx q[0];
rz(0.32904539) q[0];
rz(-1.344205) q[1];
sx q[1];
rz(-2.3842594) q[1];
sx q[1];
rz(3.0175041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24606516) q[0];
sx q[0];
rz(-2.82282) q[0];
sx q[0];
rz(0.84455873) q[0];
rz(2.7001405) q[2];
sx q[2];
rz(-2.8270671) q[2];
sx q[2];
rz(-2.0528094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.263843) q[1];
sx q[1];
rz(-1.2141435) q[1];
sx q[1];
rz(2.043173) q[1];
rz(-0.20601087) q[3];
sx q[3];
rz(-2.2643746) q[3];
sx q[3];
rz(-1.8172952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2825534) q[2];
sx q[2];
rz(-0.48290792) q[2];
sx q[2];
rz(-0.60749751) q[2];
rz(2.2533158) q[3];
sx q[3];
rz(-1.6162623) q[3];
sx q[3];
rz(1.4265149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68179503) q[0];
sx q[0];
rz(-1.2920222) q[0];
sx q[0];
rz(-0.23455308) q[0];
rz(1.0598496) q[1];
sx q[1];
rz(-1.1666965) q[1];
sx q[1];
rz(0.92811981) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1477032) q[0];
sx q[0];
rz(-0.77489955) q[0];
sx q[0];
rz(-2.5323917) q[0];
x q[1];
rz(1.9869204) q[2];
sx q[2];
rz(-2.0095996) q[2];
sx q[2];
rz(3.0896387) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8923087) q[1];
sx q[1];
rz(-1.4362446) q[1];
sx q[1];
rz(-2.8110912) q[1];
rz(0.15367963) q[3];
sx q[3];
rz(-1.9369594) q[3];
sx q[3];
rz(1.8813713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5633391) q[2];
sx q[2];
rz(-1.7208865) q[2];
sx q[2];
rz(1.2467747) q[2];
rz(-2.9747544) q[3];
sx q[3];
rz(-2.2127547) q[3];
sx q[3];
rz(1.6125352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.700915) q[0];
sx q[0];
rz(-0.066078521) q[0];
sx q[0];
rz(2.3833158) q[0];
rz(1.5757163) q[1];
sx q[1];
rz(-2.5570452) q[1];
sx q[1];
rz(-0.87055269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2038531) q[0];
sx q[0];
rz(-2.014262) q[0];
sx q[0];
rz(0.70758836) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72239082) q[2];
sx q[2];
rz(-0.71522994) q[2];
sx q[2];
rz(-1.720429) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16337559) q[1];
sx q[1];
rz(-1.723147) q[1];
sx q[1];
rz(-1.9780897) q[1];
rz(-2.5136886) q[3];
sx q[3];
rz(-1.5122754) q[3];
sx q[3];
rz(-2.1128775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7190651) q[2];
sx q[2];
rz(-1.0079404) q[2];
sx q[2];
rz(-2.3818805) q[2];
rz(0.64940137) q[3];
sx q[3];
rz(-1.6067182) q[3];
sx q[3];
rz(-1.5627129) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11288697) q[0];
sx q[0];
rz(-0.6539456) q[0];
sx q[0];
rz(1.9586067) q[0];
rz(0.68823367) q[1];
sx q[1];
rz(-1.3166683) q[1];
sx q[1];
rz(0.36929718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6128176) q[0];
sx q[0];
rz(-1.5982096) q[0];
sx q[0];
rz(-0.66365906) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8393272) q[2];
sx q[2];
rz(-2.4068953) q[2];
sx q[2];
rz(2.1718028) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10864133) q[1];
sx q[1];
rz(-0.88283112) q[1];
sx q[1];
rz(0.23984075) q[1];
rz(-pi) q[2];
rz(1.9019674) q[3];
sx q[3];
rz(-1.8871565) q[3];
sx q[3];
rz(1.6744194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1987622) q[2];
sx q[2];
rz(-1.3254712) q[2];
sx q[2];
rz(0.99786264) q[2];
rz(-2.5241847) q[3];
sx q[3];
rz(-2.182775) q[3];
sx q[3];
rz(0.82120419) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038641039) q[0];
sx q[0];
rz(-1.2741673) q[0];
sx q[0];
rz(-1.2432903) q[0];
rz(0.78168166) q[1];
sx q[1];
rz(-1.8676753) q[1];
sx q[1];
rz(-0.26085687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7383682) q[0];
sx q[0];
rz(-2.0595831) q[0];
sx q[0];
rz(1.8098149) q[0];
rz(-pi) q[1];
rz(-1.1356372) q[2];
sx q[2];
rz(-1.9138304) q[2];
sx q[2];
rz(-1.940286) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32567027) q[1];
sx q[1];
rz(-1.7223893) q[1];
sx q[1];
rz(-1.1071015) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2844159) q[3];
sx q[3];
rz(-2.4900511) q[3];
sx q[3];
rz(2.1809354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47238749) q[2];
sx q[2];
rz(-2.278625) q[2];
sx q[2];
rz(-0.48259398) q[2];
rz(-3.0274296) q[3];
sx q[3];
rz(-2.0518905) q[3];
sx q[3];
rz(0.94858661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.362185) q[0];
sx q[0];
rz(-0.055883378) q[0];
sx q[0];
rz(2.8756397) q[0];
rz(0.42568046) q[1];
sx q[1];
rz(-1.2454147) q[1];
sx q[1];
rz(-0.46806213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6702995) q[0];
sx q[0];
rz(-2.6001996) q[0];
sx q[0];
rz(1.584021) q[0];
rz(-1.9114248) q[2];
sx q[2];
rz(-1.4832895) q[2];
sx q[2];
rz(2.3208095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7167042) q[1];
sx q[1];
rz(-2.0603097) q[1];
sx q[1];
rz(-2.8498226) q[1];
rz(-pi) q[2];
rz(-2.0909833) q[3];
sx q[3];
rz(-1.597763) q[3];
sx q[3];
rz(0.20668465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.88428664) q[2];
sx q[2];
rz(-1.2904737) q[2];
sx q[2];
rz(-0.5536983) q[2];
rz(-0.92343679) q[3];
sx q[3];
rz(-2.2348576) q[3];
sx q[3];
rz(-0.060801774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629352) q[0];
sx q[0];
rz(-2.904992) q[0];
sx q[0];
rz(2.4635354) q[0];
rz(2.1642115) q[1];
sx q[1];
rz(-1.1481608) q[1];
sx q[1];
rz(1.703702) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2861833) q[0];
sx q[0];
rz(-2.0269505) q[0];
sx q[0];
rz(2.9265704) q[0];
x q[1];
rz(1.5680268) q[2];
sx q[2];
rz(-1.0068147) q[2];
sx q[2];
rz(1.1344879) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9919093) q[1];
sx q[1];
rz(-1.339773) q[1];
sx q[1];
rz(2.8808589) q[1];
rz(-pi) q[2];
rz(2.3947761) q[3];
sx q[3];
rz(-1.9821315) q[3];
sx q[3];
rz(0.29838994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.73072726) q[2];
sx q[2];
rz(-0.56034708) q[2];
sx q[2];
rz(-0.66696683) q[2];
rz(3.0242331) q[3];
sx q[3];
rz(-1.8662235) q[3];
sx q[3];
rz(1.3807152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.9475107) q[0];
sx q[0];
rz(-1.5503333) q[0];
sx q[0];
rz(-2.7981753) q[0];
rz(-1.4944448) q[1];
sx q[1];
rz(-2.1518555) q[1];
sx q[1];
rz(0.48430482) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0065932) q[0];
sx q[0];
rz(-1.8350701) q[0];
sx q[0];
rz(-2.5359383) q[0];
x q[1];
rz(1.1555919) q[2];
sx q[2];
rz(-1.9671429) q[2];
sx q[2];
rz(1.1434778) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2769988) q[1];
sx q[1];
rz(-0.62886695) q[1];
sx q[1];
rz(0.38284812) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4196641) q[3];
sx q[3];
rz(-0.32708229) q[3];
sx q[3];
rz(2.7308635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9453498) q[2];
sx q[2];
rz(-1.1002772) q[2];
sx q[2];
rz(2.2958882) q[2];
rz(1.9218933) q[3];
sx q[3];
rz(-1.2518576) q[3];
sx q[3];
rz(-2.9367101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4591111) q[0];
sx q[0];
rz(-2.8963608) q[0];
sx q[0];
rz(2.5526175) q[0];
rz(2.466195) q[1];
sx q[1];
rz(-2.1928936) q[1];
sx q[1];
rz(-1.677547) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16892085) q[0];
sx q[0];
rz(-1.1612478) q[0];
sx q[0];
rz(0.99486339) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.079564261) q[2];
sx q[2];
rz(-0.83919385) q[2];
sx q[2];
rz(-2.6688719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3879536) q[1];
sx q[1];
rz(-1.6170701) q[1];
sx q[1];
rz(-0.58677499) q[1];
x q[2];
rz(2.6028529) q[3];
sx q[3];
rz(-2.3711088) q[3];
sx q[3];
rz(-2.89794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4097164) q[2];
sx q[2];
rz(-1.3940553) q[2];
sx q[2];
rz(-2.0173006) q[2];
rz(-2.2935947) q[3];
sx q[3];
rz(-2.336899) q[3];
sx q[3];
rz(-3.1112352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59415862) q[0];
sx q[0];
rz(-1.1347329) q[0];
sx q[0];
rz(-1.5190079) q[0];
rz(-0.92319725) q[1];
sx q[1];
rz(-2.1122439) q[1];
sx q[1];
rz(-1.5069638) q[1];
rz(2.6106264) q[2];
sx q[2];
rz(-2.785977) q[2];
sx q[2];
rz(-1.2422864) q[2];
rz(1.2673479) q[3];
sx q[3];
rz(-1.787723) q[3];
sx q[3];
rz(2.1026305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
