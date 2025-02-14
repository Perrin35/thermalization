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
rz(-0.75495523) q[0];
sx q[0];
rz(-2.39019) q[0];
sx q[0];
rz(1.0487392) q[0];
rz(0.41930786) q[1];
sx q[1];
rz(-1.3860621) q[1];
sx q[1];
rz(0.11828932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8784507) q[0];
sx q[0];
rz(-1.0174482) q[0];
sx q[0];
rz(2.8216777) q[0];
rz(-1.5544791) q[2];
sx q[2];
rz(-1.6217578) q[2];
sx q[2];
rz(1.9555443) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1227545) q[1];
sx q[1];
rz(-1.600482) q[1];
sx q[1];
rz(1.5537986) q[1];
rz(-pi) q[2];
rz(-2.8382073) q[3];
sx q[3];
rz(-1.4974563) q[3];
sx q[3];
rz(-0.23305063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92363858) q[2];
sx q[2];
rz(-3.1380234) q[2];
sx q[2];
rz(0.16036073) q[2];
rz(1.977836) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639975) q[0];
sx q[0];
rz(-1.4871335) q[0];
sx q[0];
rz(0.98597041) q[0];
rz(1.5522955) q[1];
sx q[1];
rz(-2.8542216) q[1];
sx q[1];
rz(1.5602962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6628806) q[0];
sx q[0];
rz(-1.1887738) q[0];
sx q[0];
rz(-2.0860614) q[0];
x q[1];
rz(2.5471117) q[2];
sx q[2];
rz(-1.5533326) q[2];
sx q[2];
rz(1.6169777) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.040890732) q[1];
sx q[1];
rz(-1.9305848) q[1];
sx q[1];
rz(-1.3967819) q[1];
rz(-pi) q[2];
rz(-1.738882) q[3];
sx q[3];
rz(-2.5017284) q[3];
sx q[3];
rz(2.9652023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61957773) q[2];
sx q[2];
rz(-1.8176983) q[2];
sx q[2];
rz(2.1066693) q[2];
rz(-2.6400635) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(0.97626221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72306776) q[0];
sx q[0];
rz(-1.9880966) q[0];
sx q[0];
rz(-1.9368517) q[0];
rz(-1.0459666) q[1];
sx q[1];
rz(-3.0515262) q[1];
sx q[1];
rz(0.14030309) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2818272) q[0];
sx q[0];
rz(-0.66376309) q[0];
sx q[0];
rz(-1.1361994) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42201747) q[2];
sx q[2];
rz(-0.8895424) q[2];
sx q[2];
rz(-1.2764507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0187776) q[1];
sx q[1];
rz(-2.5053456) q[1];
sx q[1];
rz(2.754209) q[1];
x q[2];
rz(1.9123593) q[3];
sx q[3];
rz(-1.7031125) q[3];
sx q[3];
rz(1.7225456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3673765) q[2];
sx q[2];
rz(-0.99952951) q[2];
sx q[2];
rz(2.9191169) q[2];
rz(-0.065464822) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(-2.2953667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881994) q[0];
sx q[0];
rz(-3.1173752) q[0];
sx q[0];
rz(-0.54704332) q[0];
rz(-0.25829265) q[1];
sx q[1];
rz(-3.1196085) q[1];
sx q[1];
rz(-2.8000854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40453005) q[0];
sx q[0];
rz(-1.5642691) q[0];
sx q[0];
rz(-3.1402236) q[0];
rz(-pi) q[1];
rz(1.1522305) q[2];
sx q[2];
rz(-1.8822877) q[2];
sx q[2];
rz(-0.77889393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66879067) q[1];
sx q[1];
rz(-0.83757949) q[1];
sx q[1];
rz(-1.1146783) q[1];
x q[2];
rz(-0.82587256) q[3];
sx q[3];
rz(-2.4950728) q[3];
sx q[3];
rz(-0.41057983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7121938) q[2];
sx q[2];
rz(-1.2960351) q[2];
sx q[2];
rz(-2.1923501) q[2];
rz(2.3927169) q[3];
sx q[3];
rz(-1.2810992) q[3];
sx q[3];
rz(-0.062189814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.6147989) q[0];
sx q[0];
rz(-3.1068046) q[0];
sx q[0];
rz(-1.5507966) q[0];
rz(-1.7855478) q[1];
sx q[1];
rz(-3.1372012) q[1];
sx q[1];
rz(-0.063025085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9585091) q[0];
sx q[0];
rz(-1.5048072) q[0];
sx q[0];
rz(-1.5404758) q[0];
rz(0.75205113) q[2];
sx q[2];
rz(-1.1630485) q[2];
sx q[2];
rz(-2.4501981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5668727) q[1];
sx q[1];
rz(-1.5441455) q[1];
sx q[1];
rz(-1.779056) q[1];
rz(-pi) q[2];
rz(1.0980561) q[3];
sx q[3];
rz(-2.2946828) q[3];
sx q[3];
rz(-0.19886097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65681347) q[2];
sx q[2];
rz(-1.8317089) q[2];
sx q[2];
rz(2.4978034) q[2];
rz(2.2667609) q[3];
sx q[3];
rz(-0.30431408) q[3];
sx q[3];
rz(-2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0952045) q[0];
sx q[0];
rz(-0.055618532) q[0];
sx q[0];
rz(2.7860506) q[0];
rz(0.19861673) q[1];
sx q[1];
rz(-0.0067409975) q[1];
sx q[1];
rz(2.9933062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3902238) q[0];
sx q[0];
rz(-1.7538188) q[0];
sx q[0];
rz(1.5682632) q[0];
x q[1];
rz(2.9695187) q[2];
sx q[2];
rz(-0.41671696) q[2];
sx q[2];
rz(-0.91815776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4365789) q[1];
sx q[1];
rz(-1.8039898) q[1];
sx q[1];
rz(1.4326653) q[1];
rz(-pi) q[2];
rz(-1.8174174) q[3];
sx q[3];
rz(-2.518836) q[3];
sx q[3];
rz(2.0758219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6716914) q[2];
sx q[2];
rz(-2.9000059) q[2];
sx q[2];
rz(0.12413231) q[2];
rz(2.5668674) q[3];
sx q[3];
rz(-0.14437965) q[3];
sx q[3];
rz(-2.9706484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8827051) q[0];
sx q[0];
rz(-0.12508617) q[0];
sx q[0];
rz(0.74147725) q[0];
rz(0.28400907) q[1];
sx q[1];
rz(-0.0037071204) q[1];
sx q[1];
rz(0.31518087) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18703546) q[0];
sx q[0];
rz(-1.5976359) q[0];
sx q[0];
rz(-1.5052273) q[0];
x q[1];
rz(2.0186606) q[2];
sx q[2];
rz(-1.383198) q[2];
sx q[2];
rz(-2.5459187) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6446723) q[1];
sx q[1];
rz(-2.5240891) q[1];
sx q[1];
rz(-0.27568494) q[1];
x q[2];
rz(0.02353286) q[3];
sx q[3];
rz(-1.3472024) q[3];
sx q[3];
rz(1.6361039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2726941) q[2];
sx q[2];
rz(-1.0947451) q[2];
sx q[2];
rz(0.72186738) q[2];
rz(0.38247153) q[3];
sx q[3];
rz(-1.9964652) q[3];
sx q[3];
rz(-2.0731879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5751936) q[0];
sx q[0];
rz(-0.02481758) q[0];
sx q[0];
rz(-1.5665293) q[0];
rz(2.9381835) q[1];
sx q[1];
rz(-1.8433488) q[1];
sx q[1];
rz(0.64483109) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0087850182) q[0];
sx q[0];
rz(-1.5791348) q[0];
sx q[0];
rz(1.3632644) q[0];
rz(-0.036083607) q[2];
sx q[2];
rz(-2.5390194) q[2];
sx q[2];
rz(2.7831603) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37721021) q[1];
sx q[1];
rz(-1.4971716) q[1];
sx q[1];
rz(-1.0318569) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6504114) q[3];
sx q[3];
rz(-1.1648263) q[3];
sx q[3];
rz(0.58455672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5975534) q[2];
sx q[2];
rz(-0.35352239) q[2];
sx q[2];
rz(0.37975797) q[2];
rz(-2.0960268) q[3];
sx q[3];
rz(-1.9087722) q[3];
sx q[3];
rz(-1.9426965) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7757292) q[0];
sx q[0];
rz(-0.033670306) q[0];
sx q[0];
rz(-1.7849543) q[0];
rz(2.7011073) q[1];
sx q[1];
rz(-2.0511274) q[1];
sx q[1];
rz(0.7007362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7788575) q[0];
sx q[0];
rz(-2.252251) q[0];
sx q[0];
rz(-0.89618857) q[0];
rz(-1.8092625) q[2];
sx q[2];
rz(-2.4372413) q[2];
sx q[2];
rz(2.4865502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0807018) q[1];
sx q[1];
rz(-1.3403088) q[1];
sx q[1];
rz(-2.4031765) q[1];
rz(-pi) q[2];
rz(0.0088720284) q[3];
sx q[3];
rz(-2.7089467) q[3];
sx q[3];
rz(0.052730058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35357722) q[2];
sx q[2];
rz(-0.37166301) q[2];
sx q[2];
rz(-1.2865944) q[2];
rz(-0.50518099) q[3];
sx q[3];
rz(-2.6954539) q[3];
sx q[3];
rz(1.2222458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.5196359) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(1.5420445) q[0];
rz(2.385335) q[1];
sx q[1];
rz(-3.1344963) q[1];
sx q[1];
rz(-0.33682987) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3339506) q[0];
sx q[0];
rz(-1.5686638) q[0];
sx q[0];
rz(2.8239408) q[0];
rz(-pi) q[1];
rz(-0.9438528) q[2];
sx q[2];
rz(-0.34239951) q[2];
sx q[2];
rz(-1.8998208) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.017180786) q[1];
sx q[1];
rz(-1.4596816) q[1];
sx q[1];
rz(3.0558056) q[1];
rz(-pi) q[2];
rz(0.95864166) q[3];
sx q[3];
rz(-2.4869182) q[3];
sx q[3];
rz(-2.502113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51356641) q[2];
sx q[2];
rz(-0.94945532) q[2];
sx q[2];
rz(-0.27964082) q[2];
rz(-2.5102992) q[3];
sx q[3];
rz(-2.203233) q[3];
sx q[3];
rz(-0.62425557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8414128) q[0];
sx q[0];
rz(-1.5914088) q[0];
sx q[0];
rz(1.7803022) q[0];
rz(-0.77371669) q[1];
sx q[1];
rz(-2.5061889) q[1];
sx q[1];
rz(-2.9255964) q[1];
rz(0.85971268) q[2];
sx q[2];
rz(-2.0794686) q[2];
sx q[2];
rz(2.5735264) q[2];
rz(-1.5935043) q[3];
sx q[3];
rz(-1.9469713) q[3];
sx q[3];
rz(3.136829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
