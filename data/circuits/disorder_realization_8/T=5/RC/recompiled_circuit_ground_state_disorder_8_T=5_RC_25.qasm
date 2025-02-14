OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.79787624) q[0];
sx q[0];
rz(-0.98348445) q[0];
sx q[0];
rz(2.9498192) q[0];
rz(0.51044381) q[1];
sx q[1];
rz(4.508701) q[1];
sx q[1];
rz(8.0291168) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5001427) q[0];
sx q[0];
rz(-0.035273835) q[0];
sx q[0];
rz(-2.7462237) q[0];
rz(-2.4641855) q[2];
sx q[2];
rz(-2.8500994) q[2];
sx q[2];
rz(-0.10805932) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6916445) q[1];
sx q[1];
rz(-2.3385923) q[1];
sx q[1];
rz(1.8464441) q[1];
rz(-pi) q[2];
rz(-1.1825425) q[3];
sx q[3];
rz(-2.7388529) q[3];
sx q[3];
rz(-2.5975063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.50769606) q[2];
sx q[2];
rz(-2.8140929) q[2];
sx q[2];
rz(-0.92602777) q[2];
rz(-1.6294468) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(-1.1431471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73285) q[0];
sx q[0];
rz(-0.73246211) q[0];
sx q[0];
rz(2.9090885) q[0];
rz(-2.0022424) q[1];
sx q[1];
rz(-1.573223) q[1];
sx q[1];
rz(2.2154636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0588603) q[0];
sx q[0];
rz(-1.345536) q[0];
sx q[0];
rz(-2.9316069) q[0];
rz(-pi) q[1];
rz(-2.813999) q[2];
sx q[2];
rz(-1.4783876) q[2];
sx q[2];
rz(1.7756697) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7826816) q[1];
sx q[1];
rz(-2.5645683) q[1];
sx q[1];
rz(-0.38459528) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1180463) q[3];
sx q[3];
rz(-1.0501852) q[3];
sx q[3];
rz(1.758616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53655475) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(-0.79263318) q[2];
rz(0.41147453) q[3];
sx q[3];
rz(-2.1992407) q[3];
sx q[3];
rz(-2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8234392) q[0];
sx q[0];
rz(-2.1397488) q[0];
sx q[0];
rz(-2.8807628) q[0];
rz(-2.8584495) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(1.0556489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94819647) q[0];
sx q[0];
rz(-2.4760111) q[0];
sx q[0];
rz(-1.2917305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63313578) q[2];
sx q[2];
rz(-1.5738259) q[2];
sx q[2];
rz(-2.5315447) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4312517) q[1];
sx q[1];
rz(-1.6762937) q[1];
sx q[1];
rz(-0.72781095) q[1];
rz(-pi) q[2];
rz(1.0751455) q[3];
sx q[3];
rz(-1.6248253) q[3];
sx q[3];
rz(3.0040405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5259214) q[2];
sx q[2];
rz(-1.0893704) q[2];
sx q[2];
rz(1.2582568) q[2];
rz(-1.0387756) q[3];
sx q[3];
rz(-2.153331) q[3];
sx q[3];
rz(-1.4163777) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0106169) q[0];
sx q[0];
rz(-1.6058141) q[0];
sx q[0];
rz(-3.0220939) q[0];
rz(-0.15826982) q[1];
sx q[1];
rz(-0.4522849) q[1];
sx q[1];
rz(-1.7012885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8058469) q[0];
sx q[0];
rz(-0.81253883) q[0];
sx q[0];
rz(-1.0643908) q[0];
rz(1.2358627) q[2];
sx q[2];
rz(-1.6455629) q[2];
sx q[2];
rz(2.3521545) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56934565) q[1];
sx q[1];
rz(-1.3075446) q[1];
sx q[1];
rz(-2.6514228) q[1];
rz(0.3500895) q[3];
sx q[3];
rz(-2.0357657) q[3];
sx q[3];
rz(1.657682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4629472) q[2];
sx q[2];
rz(-3.0061649) q[2];
sx q[2];
rz(0.43369183) q[2];
rz(-0.67484754) q[3];
sx q[3];
rz(-1.0899455) q[3];
sx q[3];
rz(-1.1564144) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40815142) q[0];
sx q[0];
rz(-2.0858522) q[0];
sx q[0];
rz(2.5422886) q[0];
rz(0.76404461) q[1];
sx q[1];
rz(-2.1115477) q[1];
sx q[1];
rz(-2.5291671) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0611506) q[0];
sx q[0];
rz(-2.168405) q[0];
sx q[0];
rz(-2.963128) q[0];
rz(-3.1407498) q[2];
sx q[2];
rz(-0.71686059) q[2];
sx q[2];
rz(-1.5886024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45578411) q[1];
sx q[1];
rz(-2.1791885) q[1];
sx q[1];
rz(1.366282) q[1];
x q[2];
rz(2.3633943) q[3];
sx q[3];
rz(-0.68420568) q[3];
sx q[3];
rz(-2.9915031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4446438) q[2];
sx q[2];
rz(-0.79978839) q[2];
sx q[2];
rz(-2.6714228) q[2];
rz(-2.7759975) q[3];
sx q[3];
rz(-1.1919034) q[3];
sx q[3];
rz(2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81569165) q[0];
sx q[0];
rz(-0.79623643) q[0];
sx q[0];
rz(0.66194397) q[0];
rz(0.081261948) q[1];
sx q[1];
rz(-0.47835246) q[1];
sx q[1];
rz(3.001396) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79665444) q[0];
sx q[0];
rz(-1.7895763) q[0];
sx q[0];
rz(-0.47056912) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1344994) q[2];
sx q[2];
rz(-2.4738174) q[2];
sx q[2];
rz(0.89950022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0602136) q[1];
sx q[1];
rz(-1.5790931) q[1];
sx q[1];
rz(-1.0992728) q[1];
x q[2];
rz(1.3883038) q[3];
sx q[3];
rz(-1.2267707) q[3];
sx q[3];
rz(-1.1870015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36416546) q[2];
sx q[2];
rz(-1.825793) q[2];
sx q[2];
rz(2.8167456) q[2];
rz(-2.9835564) q[3];
sx q[3];
rz(-0.41301781) q[3];
sx q[3];
rz(-2.1260927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3196816) q[0];
sx q[0];
rz(-3.0653636) q[0];
sx q[0];
rz(2.2538189) q[0];
rz(-2.7582788) q[1];
sx q[1];
rz(-2.2593468) q[1];
sx q[1];
rz(-2.9820138) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5349284) q[0];
sx q[0];
rz(-1.3114358) q[0];
sx q[0];
rz(1.242437) q[0];
x q[1];
rz(0.55241199) q[2];
sx q[2];
rz(-1.0670976) q[2];
sx q[2];
rz(-2.9744862) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23160203) q[1];
sx q[1];
rz(-1.9456777) q[1];
sx q[1];
rz(1.433062) q[1];
x q[2];
rz(0.28954808) q[3];
sx q[3];
rz(-1.411045) q[3];
sx q[3];
rz(-1.8600132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5923656) q[2];
sx q[2];
rz(-1.9840019) q[2];
sx q[2];
rz(-1.7249736) q[2];
rz(1.8727411) q[3];
sx q[3];
rz(-2.3609991) q[3];
sx q[3];
rz(-2.1835073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6855327) q[0];
sx q[0];
rz(-1.1153509) q[0];
sx q[0];
rz(0.90150315) q[0];
rz(-0.69425663) q[1];
sx q[1];
rz(-1.4638823) q[1];
sx q[1];
rz(-1.136397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89595786) q[0];
sx q[0];
rz(-0.16366009) q[0];
sx q[0];
rz(1.3748522) q[0];
rz(-pi) q[1];
rz(-0.56635277) q[2];
sx q[2];
rz(-2.512062) q[2];
sx q[2];
rz(0.16816631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7982805) q[1];
sx q[1];
rz(-1.5250051) q[1];
sx q[1];
rz(-0.87516038) q[1];
rz(-pi) q[2];
rz(2.091892) q[3];
sx q[3];
rz(-0.9113833) q[3];
sx q[3];
rz(-1.570861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14094341) q[2];
sx q[2];
rz(-1.2691754) q[2];
sx q[2];
rz(-2.3285274) q[2];
rz(0.55417577) q[3];
sx q[3];
rz(-1.7289836) q[3];
sx q[3];
rz(0.36177844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.5799705) q[0];
sx q[0];
rz(-1.1142092) q[0];
sx q[0];
rz(2.9680874) q[0];
rz(2.7165727) q[1];
sx q[1];
rz(-1.5360906) q[1];
sx q[1];
rz(0.82957155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3542784) q[0];
sx q[0];
rz(-1.4256546) q[0];
sx q[0];
rz(0.73696846) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.095171) q[2];
sx q[2];
rz(-1.5299462) q[2];
sx q[2];
rz(-2.61039) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8977938) q[1];
sx q[1];
rz(-1.3659371) q[1];
sx q[1];
rz(-2.5216991) q[1];
rz(-pi) q[2];
rz(-0.18316571) q[3];
sx q[3];
rz(-0.34160994) q[3];
sx q[3];
rz(-0.46942018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7912264) q[2];
sx q[2];
rz(-1.841505) q[2];
sx q[2];
rz(-1.4607956) q[2];
rz(-0.99460498) q[3];
sx q[3];
rz(-0.9959144) q[3];
sx q[3];
rz(-2.2554956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17722546) q[0];
sx q[0];
rz(-0.57562861) q[0];
sx q[0];
rz(-3.0395569) q[0];
rz(0.61965865) q[1];
sx q[1];
rz(-1.4980059) q[1];
sx q[1];
rz(-0.59725753) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16011691) q[0];
sx q[0];
rz(-0.61321027) q[0];
sx q[0];
rz(-0.23399467) q[0];
rz(-pi) q[1];
rz(0.37924699) q[2];
sx q[2];
rz(-2.7617117) q[2];
sx q[2];
rz(2.3178562) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1015375) q[1];
sx q[1];
rz(-1.005508) q[1];
sx q[1];
rz(0.043626309) q[1];
rz(-2.6701376) q[3];
sx q[3];
rz(-1.5999891) q[3];
sx q[3];
rz(-1.7144793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9336046) q[2];
sx q[2];
rz(-0.36078578) q[2];
sx q[2];
rz(2.2130845) q[2];
rz(1.1601296) q[3];
sx q[3];
rz(-1.4312294) q[3];
sx q[3];
rz(0.76735705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3930897) q[0];
sx q[0];
rz(-2.4414283) q[0];
sx q[0];
rz(0.23165942) q[0];
rz(-1.1275935) q[1];
sx q[1];
rz(-0.53378202) q[1];
sx q[1];
rz(3.1376874) q[1];
rz(1.4814346) q[2];
sx q[2];
rz(-0.67831525) q[2];
sx q[2];
rz(-0.44616551) q[2];
rz(1.9740022) q[3];
sx q[3];
rz(-0.73763631) q[3];
sx q[3];
rz(0.49745001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
