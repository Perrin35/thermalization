OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(2.0895045) q[0];
sx q[0];
rz(10.917561) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(5.5881349) q[1];
sx q[1];
rz(10.499428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.933691) q[0];
sx q[0];
rz(-1.5760734) q[0];
sx q[0];
rz(3.123377) q[0];
rz(-pi) q[1];
rz(-1.7311814) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(-0.10057848) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97779146) q[1];
sx q[1];
rz(-1.4492387) q[1];
sx q[1];
rz(1.4896643) q[1];
x q[2];
rz(-0.69016506) q[3];
sx q[3];
rz(-2.4881722) q[3];
sx q[3];
rz(1.1284459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82912123) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.5134229) q[2];
rz(-1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-0.2127969) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56869498) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(-2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.5011903) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6340856) q[0];
sx q[0];
rz(-1.5956143) q[0];
sx q[0];
rz(-3.0188942) q[0];
rz(-pi) q[1];
rz(-2.8295849) q[2];
sx q[2];
rz(-2.0199676) q[2];
sx q[2];
rz(2.1829407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7306819) q[1];
sx q[1];
rz(-0.69697471) q[1];
sx q[1];
rz(0.26941401) q[1];
rz(-pi) q[2];
rz(-0.0060175671) q[3];
sx q[3];
rz(-1.9618417) q[3];
sx q[3];
rz(1.0139549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.286065) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(-1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062155) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(1.9619933) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(-1.5037781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5119748) q[0];
sx q[0];
rz(-3.0324728) q[0];
sx q[0];
rz(2.2114121) q[0];
x q[1];
rz(0.47567993) q[2];
sx q[2];
rz(-2.1026346) q[2];
sx q[2];
rz(1.4059517) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40239247) q[1];
sx q[1];
rz(-1.270073) q[1];
sx q[1];
rz(-1.5261569) q[1];
rz(-pi) q[2];
rz(0.13111968) q[3];
sx q[3];
rz(-0.91851202) q[3];
sx q[3];
rz(-0.91344792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(-0.25742325) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1733615) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(1.051735) q[0];
rz(2.5412718) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.9925041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79608941) q[0];
sx q[0];
rz(-2.8755099) q[0];
sx q[0];
rz(2.7059113) q[0];
rz(-pi) q[1];
rz(-0.97781424) q[2];
sx q[2];
rz(-1.2056418) q[2];
sx q[2];
rz(-2.1821373) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1059873) q[1];
sx q[1];
rz(-0.78250256) q[1];
sx q[1];
rz(1.1483907) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37255128) q[3];
sx q[3];
rz(-2.3636892) q[3];
sx q[3];
rz(2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(0.66876283) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4887061) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(1.0523798) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(0.88422424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54908862) q[0];
sx q[0];
rz(-1.3929402) q[0];
sx q[0];
rz(-1.7643719) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82489478) q[2];
sx q[2];
rz(-2.5755304) q[2];
sx q[2];
rz(2.538946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4864038) q[1];
sx q[1];
rz(-0.15895325) q[1];
sx q[1];
rz(2.7093191) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36966427) q[3];
sx q[3];
rz(-2.0737994) q[3];
sx q[3];
rz(-0.30077416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23652442) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(1.0926251) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
rz(1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.3177692) q[0];
rz(-0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7965308) q[0];
sx q[0];
rz(-3.025265) q[0];
sx q[0];
rz(1.2141114) q[0];
rz(-1.5029807) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(0.80578795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8816996) q[1];
sx q[1];
rz(-1.1105781) q[1];
sx q[1];
rz(-2.0395181) q[1];
rz(-pi) q[2];
x q[2];
rz(-3*pi/11) q[3];
sx q[3];
rz(-2.487605) q[3];
sx q[3];
rz(1.6884782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51450729) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(-0.43760854) q[2];
rz(-0.058241025) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(1.4490022) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.4356027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5098269) q[0];
sx q[0];
rz(-2.0938211) q[0];
sx q[0];
rz(2.564389) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1232386) q[2];
sx q[2];
rz(-1.5905426) q[2];
sx q[2];
rz(0.28723082) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0457397) q[1];
sx q[1];
rz(-0.47912712) q[1];
sx q[1];
rz(-2.2687001) q[1];
x q[2];
rz(-0.42922677) q[3];
sx q[3];
rz(-2.0373404) q[3];
sx q[3];
rz(0.45688094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-0.1329578) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4505969) q[0];
sx q[0];
rz(-0.15582514) q[0];
sx q[0];
rz(2.0777006) q[0];
x q[1];
rz(-1.442996) q[2];
sx q[2];
rz(-1.5926077) q[2];
sx q[2];
rz(-0.70805659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8088088) q[1];
sx q[1];
rz(-0.53408264) q[1];
sx q[1];
rz(-3.0384455) q[1];
x q[2];
rz(-2.7216464) q[3];
sx q[3];
rz(-1.4010324) q[3];
sx q[3];
rz(1.1743634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0662971) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(-0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(2.8310006) q[0];
rz(2.8839135) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(2.2176567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7565219) q[0];
sx q[0];
rz(-2.0002881) q[0];
sx q[0];
rz(-2.4247652) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1191101) q[2];
sx q[2];
rz(-2.6528931) q[2];
sx q[2];
rz(1.894001) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3194771) q[1];
sx q[1];
rz(-1.176782) q[1];
sx q[1];
rz(1.0072717) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68848916) q[3];
sx q[3];
rz(-0.8430891) q[3];
sx q[3];
rz(-1.9635995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6442287) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(-1.3941992) q[2];
rz(-1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7681463) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-0.038169233) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.569012) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(-1.6174181) q[0];
x q[1];
rz(-1.795904) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(-3.0181146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.888615) q[1];
sx q[1];
rz(-1.9441609) q[1];
sx q[1];
rz(3.0967767) q[1];
rz(0.90922728) q[3];
sx q[3];
rz(-1.070676) q[3];
sx q[3];
rz(1.3487032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.9434628) q[2];
rz(1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643628) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(-1.6434796) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(1.5126363) q[2];
sx q[2];
rz(-1.6806316) q[2];
sx q[2];
rz(-3.1113839) q[2];
rz(1.2002767) q[3];
sx q[3];
rz(-2.4051718) q[3];
sx q[3];
rz(-2.6889599) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
