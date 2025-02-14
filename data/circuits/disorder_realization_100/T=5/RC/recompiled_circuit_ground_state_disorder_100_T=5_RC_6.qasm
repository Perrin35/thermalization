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
rz(2.1821238) q[0];
sx q[0];
rz(9.9766599) q[0];
rz(-0.38504398) q[1];
sx q[1];
rz(4.9209891) q[1];
sx q[1];
rz(9.8434386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075567632) q[0];
sx q[0];
rz(-2.2097144) q[0];
sx q[0];
rz(0.35982168) q[0];
rz(-pi) q[1];
rz(2.6919882) q[2];
sx q[2];
rz(-1.7309963) q[2];
sx q[2];
rz(0.9148324) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6578411) q[1];
sx q[1];
rz(-2.9194909) q[1];
sx q[1];
rz(-2.7276843) q[1];
rz(-pi) q[2];
rz(1.6012264) q[3];
sx q[3];
rz(-1.5879027) q[3];
sx q[3];
rz(-2.3612983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5358413) q[2];
sx q[2];
rz(-1.8391106) q[2];
sx q[2];
rz(0.38899404) q[2];
rz(-0.89748663) q[3];
sx q[3];
rz(-0.56098452) q[3];
sx q[3];
rz(-1.2787904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9101343) q[0];
sx q[0];
rz(-0.95160216) q[0];
sx q[0];
rz(-2.8125473) q[0];
rz(1.344205) q[1];
sx q[1];
rz(-0.75733328) q[1];
sx q[1];
rz(-0.12408852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50590512) q[0];
sx q[0];
rz(-1.3342739) q[0];
sx q[0];
rz(-0.2157477) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7089073) q[2];
sx q[2];
rz(-1.8542854) q[2];
sx q[2];
rz(2.5140009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6581003) q[1];
sx q[1];
rz(-1.1303194) q[1];
sx q[1];
rz(2.7453203) q[1];
rz(-2.9355818) q[3];
sx q[3];
rz(-2.2643746) q[3];
sx q[3];
rz(-1.3242974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8590392) q[2];
sx q[2];
rz(-0.48290792) q[2];
sx q[2];
rz(-0.60749751) q[2];
rz(0.88827682) q[3];
sx q[3];
rz(-1.5253303) q[3];
sx q[3];
rz(-1.7150778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4597976) q[0];
sx q[0];
rz(-1.2920222) q[0];
sx q[0];
rz(0.23455308) q[0];
rz(-1.0598496) q[1];
sx q[1];
rz(-1.1666965) q[1];
sx q[1];
rz(2.2134728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88555777) q[0];
sx q[0];
rz(-1.9826898) q[0];
sx q[0];
rz(-0.67660634) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6675148) q[2];
sx q[2];
rz(-1.1961719) q[2];
sx q[2];
rz(1.7044391) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8923087) q[1];
sx q[1];
rz(-1.705348) q[1];
sx q[1];
rz(-2.8110912) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.987913) q[3];
sx q[3];
rz(-1.2046332) q[3];
sx q[3];
rz(-1.8813713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5633391) q[2];
sx q[2];
rz(-1.4207062) q[2];
sx q[2];
rz(1.8948179) q[2];
rz(-0.16683821) q[3];
sx q[3];
rz(-0.92883795) q[3];
sx q[3];
rz(-1.5290574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4406776) q[0];
sx q[0];
rz(-3.0755141) q[0];
sx q[0];
rz(2.3833158) q[0];
rz(-1.5658763) q[1];
sx q[1];
rz(-0.58454746) q[1];
sx q[1];
rz(0.87055269) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93773952) q[0];
sx q[0];
rz(-1.1273307) q[0];
sx q[0];
rz(-2.4340043) q[0];
rz(-pi) q[1];
rz(-0.72239082) q[2];
sx q[2];
rz(-0.71522994) q[2];
sx q[2];
rz(-1.4211637) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9782171) q[1];
sx q[1];
rz(-1.4184457) q[1];
sx q[1];
rz(1.1635029) q[1];
x q[2];
rz(1.4985256) q[3];
sx q[3];
rz(-0.94413432) q[3];
sx q[3];
rz(-0.58451239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7190651) q[2];
sx q[2];
rz(-2.1336522) q[2];
sx q[2];
rz(-0.75971216) q[2];
rz(0.64940137) q[3];
sx q[3];
rz(-1.6067182) q[3];
sx q[3];
rz(1.5788797) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11288697) q[0];
sx q[0];
rz(-2.4876471) q[0];
sx q[0];
rz(1.9586067) q[0];
rz(2.453359) q[1];
sx q[1];
rz(-1.8249244) q[1];
sx q[1];
rz(-2.7722955) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52877502) q[0];
sx q[0];
rz(-1.5982096) q[0];
sx q[0];
rz(-0.66365906) q[0];
x q[1];
rz(0.7116811) q[2];
sx q[2];
rz(-1.7717012) q[2];
sx q[2];
rz(2.7679659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0329513) q[1];
sx q[1];
rz(-0.88283112) q[1];
sx q[1];
rz(2.9017519) q[1];
rz(-pi) q[2];
rz(-0.33325382) q[3];
sx q[3];
rz(-1.8849533) q[3];
sx q[3];
rz(0.21018782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94283048) q[2];
sx q[2];
rz(-1.3254712) q[2];
sx q[2];
rz(2.14373) q[2];
rz(-2.5241847) q[3];
sx q[3];
rz(-0.95881763) q[3];
sx q[3];
rz(2.3203885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1029516) q[0];
sx q[0];
rz(-1.8674253) q[0];
sx q[0];
rz(-1.8983023) q[0];
rz(2.359911) q[1];
sx q[1];
rz(-1.8676753) q[1];
sx q[1];
rz(0.26085687) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4032245) q[0];
sx q[0];
rz(-2.0595831) q[0];
sx q[0];
rz(-1.3317778) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2736808) q[2];
sx q[2];
rz(-0.54722584) q[2];
sx q[2];
rz(0.99582129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.32567027) q[1];
sx q[1];
rz(-1.4192033) q[1];
sx q[1];
rz(2.0344911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93918899) q[3];
sx q[3];
rz(-1.3986466) q[3];
sx q[3];
rz(2.3014455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47238749) q[2];
sx q[2];
rz(-0.86296764) q[2];
sx q[2];
rz(-2.6589987) q[2];
rz(-0.11416301) q[3];
sx q[3];
rz(-1.0897021) q[3];
sx q[3];
rz(0.94858661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.362185) q[0];
sx q[0];
rz(-3.0857093) q[0];
sx q[0];
rz(2.8756397) q[0];
rz(-0.42568046) q[1];
sx q[1];
rz(-1.2454147) q[1];
sx q[1];
rz(-2.6735305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0307559) q[0];
sx q[0];
rz(-1.5776112) q[0];
sx q[0];
rz(1.0294419) q[0];
rz(-pi) q[1];
rz(-1.3139901) q[2];
sx q[2];
rz(-0.35126424) q[2];
sx q[2];
rz(-0.50826573) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99331964) q[1];
sx q[1];
rz(-2.5778505) q[1];
sx q[1];
rz(-2.0658653) q[1];
rz(-pi) q[2];
rz(-1.516582) q[3];
sx q[3];
rz(-0.52082131) q[3];
sx q[3];
rz(-1.3170742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88428664) q[2];
sx q[2];
rz(-1.2904737) q[2];
sx q[2];
rz(-2.5878944) q[2];
rz(2.2181559) q[3];
sx q[3];
rz(-2.2348576) q[3];
sx q[3];
rz(3.0807909) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629352) q[0];
sx q[0];
rz(-2.904992) q[0];
sx q[0];
rz(-0.67805725) q[0];
rz(2.1642115) q[1];
sx q[1];
rz(-1.1481608) q[1];
sx q[1];
rz(-1.4378907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7610705) q[0];
sx q[0];
rz(-1.7635404) q[0];
sx q[0];
rz(2.0362292) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5639834) q[2];
sx q[2];
rz(-1.5731369) q[2];
sx q[2];
rz(-0.43482796) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71187825) q[1];
sx q[1];
rz(-2.7949882) q[1];
sx q[1];
rz(0.73965736) q[1];
x q[2];
rz(0.74681654) q[3];
sx q[3];
rz(-1.1594611) q[3];
sx q[3];
rz(-2.8432027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73072726) q[2];
sx q[2];
rz(-2.5812456) q[2];
sx q[2];
rz(-0.66696683) q[2];
rz(3.0242331) q[3];
sx q[3];
rz(-1.8662235) q[3];
sx q[3];
rz(-1.7608775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9475107) q[0];
sx q[0];
rz(-1.5503333) q[0];
sx q[0];
rz(-2.7981753) q[0];
rz(1.4944448) q[1];
sx q[1];
rz(-2.1518555) q[1];
sx q[1];
rz(-0.48430482) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0065932) q[0];
sx q[0];
rz(-1.3065225) q[0];
sx q[0];
rz(0.60565438) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76700751) q[2];
sx q[2];
rz(-2.5756774) q[2];
sx q[2];
rz(-0.2919251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2769988) q[1];
sx q[1];
rz(-0.62886695) q[1];
sx q[1];
rz(0.38284812) q[1];
rz(-0.24933322) q[3];
sx q[3];
rz(-1.7847381) q[3];
sx q[3];
rz(0.46505022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9453498) q[2];
sx q[2];
rz(-2.0413155) q[2];
sx q[2];
rz(-2.2958882) q[2];
rz(-1.2196994) q[3];
sx q[3];
rz(-1.889735) q[3];
sx q[3];
rz(-0.20488258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.6824816) q[0];
sx q[0];
rz(-0.24523188) q[0];
sx q[0];
rz(2.5526175) q[0];
rz(2.466195) q[1];
sx q[1];
rz(-2.1928936) q[1];
sx q[1];
rz(-1.677547) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9726718) q[0];
sx q[0];
rz(-1.1612478) q[0];
sx q[0];
rz(-2.1467293) q[0];
rz(2.3039742) q[2];
sx q[2];
rz(-1.5116201) q[2];
sx q[2];
rz(2.0967332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88661609) q[1];
sx q[1];
rz(-2.5532097) q[1];
sx q[1];
rz(-3.058152) q[1];
x q[2];
rz(-2.0328224) q[3];
sx q[3];
rz(-2.2115876) q[3];
sx q[3];
rz(0.45087157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4097164) q[2];
sx q[2];
rz(-1.3940553) q[2];
sx q[2];
rz(-1.124292) q[2];
rz(2.2935947) q[3];
sx q[3];
rz(-0.8046937) q[3];
sx q[3];
rz(-3.1112352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59415862) q[0];
sx q[0];
rz(-2.0068598) q[0];
sx q[0];
rz(1.6225847) q[0];
rz(2.2183954) q[1];
sx q[1];
rz(-2.1122439) q[1];
sx q[1];
rz(-1.5069638) q[1];
rz(-2.6106264) q[2];
sx q[2];
rz(-0.35561564) q[2];
sx q[2];
rz(1.8993062) q[2];
rz(1.8742448) q[3];
sx q[3];
rz(-1.3538697) q[3];
sx q[3];
rz(-1.0389622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
