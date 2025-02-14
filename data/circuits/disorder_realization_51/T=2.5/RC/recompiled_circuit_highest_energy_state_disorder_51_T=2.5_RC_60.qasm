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
rz(3.8929953) q[0];
sx q[0];
rz(10.473517) q[0];
rz(-2.7222848) q[1];
sx q[1];
rz(4.5276548) q[1];
sx q[1];
rz(6.164896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2631419) q[0];
sx q[0];
rz(-2.1241444) q[0];
sx q[0];
rz(2.8216777) q[0];
x q[1];
rz(1.5871136) q[2];
sx q[2];
rz(-1.6217578) q[2];
sx q[2];
rz(1.9555443) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1227545) q[1];
sx q[1];
rz(-1.5411106) q[1];
sx q[1];
rz(-1.5537986) q[1];
rz(-pi) q[2];
rz(-2.8382073) q[3];
sx q[3];
rz(-1.4974563) q[3];
sx q[3];
rz(-0.23305063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2179541) q[2];
sx q[2];
rz(-3.1380234) q[2];
sx q[2];
rz(0.16036073) q[2];
rz(-1.1637566) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(-0.81419301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9775951) q[0];
sx q[0];
rz(-1.4871335) q[0];
sx q[0];
rz(-2.1556222) q[0];
rz(1.5892971) q[1];
sx q[1];
rz(-2.8542216) q[1];
sx q[1];
rz(1.5812965) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575657) q[0];
sx q[0];
rz(-2.0456893) q[0];
sx q[0];
rz(-2.7090461) q[0];
rz(-1.549717) q[2];
sx q[2];
rz(-0.97641845) q[2];
sx q[2];
rz(-3.0836058) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4680921) q[1];
sx q[1];
rz(-1.7335658) q[1];
sx q[1];
rz(2.7767608) q[1];
rz(-0.12388568) q[3];
sx q[3];
rz(-0.94138161) q[3];
sx q[3];
rz(-0.38485011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61957773) q[2];
sx q[2];
rz(-1.8176983) q[2];
sx q[2];
rz(2.1066693) q[2];
rz(0.50152913) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(0.97626221) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72306776) q[0];
sx q[0];
rz(-1.1534961) q[0];
sx q[0];
rz(1.9368517) q[0];
rz(-2.095626) q[1];
sx q[1];
rz(-0.090066411) q[1];
sx q[1];
rz(0.14030309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7810996) q[0];
sx q[0];
rz(-1.3083959) q[0];
sx q[0];
rz(0.95376063) q[0];
rz(-pi) q[1];
rz(-2.297338) q[2];
sx q[2];
rz(-1.8945969) q[2];
sx q[2];
rz(0.56991386) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5922484) q[1];
sx q[1];
rz(-2.1533433) q[1];
sx q[1];
rz(1.2986533) q[1];
rz(1.1925832) q[3];
sx q[3];
rz(-2.7762354) q[3];
sx q[3];
rz(0.20361209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3673765) q[2];
sx q[2];
rz(-0.99952951) q[2];
sx q[2];
rz(2.9191169) q[2];
rz(0.065464822) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(-0.84622598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.1881994) q[0];
sx q[0];
rz(-0.024217483) q[0];
sx q[0];
rz(0.54704332) q[0];
rz(2.8833) q[1];
sx q[1];
rz(-0.02198418) q[1];
sx q[1];
rz(-0.34150728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40453005) q[0];
sx q[0];
rz(-1.5642691) q[0];
sx q[0];
rz(0.0013690283) q[0];
rz(-pi) q[1];
rz(1.1522305) q[2];
sx q[2];
rz(-1.8822877) q[2];
sx q[2];
rz(2.3626987) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3013743) q[1];
sx q[1];
rz(-2.3010588) q[1];
sx q[1];
rz(-0.45481843) q[1];
rz(-pi) q[2];
rz(0.82587256) q[3];
sx q[3];
rz(-0.64651981) q[3];
sx q[3];
rz(2.7310128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7121938) q[2];
sx q[2];
rz(-1.2960351) q[2];
sx q[2];
rz(-0.94924259) q[2];
rz(-0.74887577) q[3];
sx q[3];
rz(-1.2810992) q[3];
sx q[3];
rz(-0.062189814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5267938) q[0];
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
rz(1.3857128) q[0];
sx q[0];
rz(-1.5405419) q[0];
sx q[0];
rz(3.0755733) q[0];
rz(0.56383987) q[2];
sx q[2];
rz(-0.83602521) q[2];
sx q[2];
rz(2.6631402) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0200796) q[1];
sx q[1];
rz(-0.20993349) q[1];
sx q[1];
rz(-1.6990183) q[1];
x q[2];
rz(1.0980561) q[3];
sx q[3];
rz(-2.2946828) q[3];
sx q[3];
rz(-0.19886097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4847792) q[2];
sx q[2];
rz(-1.8317089) q[2];
sx q[2];
rz(0.64378929) q[2];
rz(2.2667609) q[3];
sx q[3];
rz(-0.30431408) q[3];
sx q[3];
rz(-2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0463882) q[0];
sx q[0];
rz(-3.0859741) q[0];
sx q[0];
rz(-0.355542) q[0];
rz(0.19861673) q[1];
sx q[1];
rz(-3.1348517) q[1];
sx q[1];
rz(-2.9933062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3902238) q[0];
sx q[0];
rz(-1.3877739) q[0];
sx q[0];
rz(1.5682632) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9695187) q[2];
sx q[2];
rz(-0.41671696) q[2];
sx q[2];
rz(2.2234349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7050138) q[1];
sx q[1];
rz(-1.3376029) q[1];
sx q[1];
rz(-1.7089273) q[1];
rz(-1.8174174) q[3];
sx q[3];
rz(-2.518836) q[3];
sx q[3];
rz(2.0758219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6716914) q[2];
sx q[2];
rz(-0.24158676) q[2];
sx q[2];
rz(3.0174603) q[2];
rz(-0.5747253) q[3];
sx q[3];
rz(-0.14437965) q[3];
sx q[3];
rz(-2.9706484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588876) q[0];
sx q[0];
rz(-0.12508617) q[0];
sx q[0];
rz(-2.4001154) q[0];
rz(0.28400907) q[1];
sx q[1];
rz(-3.1378855) q[1];
sx q[1];
rz(2.8264118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3819987) q[0];
sx q[0];
rz(-1.505251) q[0];
sx q[0];
rz(-0.026897341) q[0];
rz(-pi) q[1];
rz(0.20756794) q[2];
sx q[2];
rz(-1.1313442) q[2];
sx q[2];
rz(1.0644827) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8410346) q[1];
sx q[1];
rz(-1.729064) q[1];
sx q[1];
rz(-2.5421418) q[1];
rz(-pi) q[2];
rz(3.1180598) q[3];
sx q[3];
rz(-1.7943903) q[3];
sx q[3];
rz(1.6361039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.86889851) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(-0.72186738) q[2];
rz(-2.7591211) q[3];
sx q[3];
rz(-1.9964652) q[3];
sx q[3];
rz(1.0684048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751936) q[0];
sx q[0];
rz(-3.1167751) q[0];
sx q[0];
rz(-1.5665293) q[0];
rz(-0.20340915) q[1];
sx q[1];
rz(-1.2982439) q[1];
sx q[1];
rz(-0.64483109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0087850182) q[0];
sx q[0];
rz(-1.5791348) q[0];
sx q[0];
rz(-1.3632644) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036083607) q[2];
sx q[2];
rz(-2.5390194) q[2];
sx q[2];
rz(-0.35843231) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7643824) q[1];
sx q[1];
rz(-1.6444211) q[1];
sx q[1];
rz(1.0318569) q[1];
x q[2];
rz(-2.6504114) q[3];
sx q[3];
rz(-1.9767663) q[3];
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
rz(-2.7618347) q[2];
rz(-2.0960268) q[3];
sx q[3];
rz(-1.2328204) q[3];
sx q[3];
rz(-1.1988962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3658635) q[0];
sx q[0];
rz(-0.033670306) q[0];
sx q[0];
rz(-1.3566383) q[0];
rz(-0.44048539) q[1];
sx q[1];
rz(-1.0904652) q[1];
sx q[1];
rz(2.4408565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3627351) q[0];
sx q[0];
rz(-0.88934169) q[0];
sx q[0];
rz(-2.2454041) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3323302) q[2];
sx q[2];
rz(-0.70435134) q[2];
sx q[2];
rz(-0.65504247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8776013) q[1];
sx q[1];
rz(-2.374568) q[1];
sx q[1];
rz(2.8061538) q[1];
rz(2.7089617) q[3];
sx q[3];
rz(-1.5670766) q[3];
sx q[3];
rz(1.631581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35357722) q[2];
sx q[2];
rz(-2.7699296) q[2];
sx q[2];
rz(1.2865944) q[2];
rz(2.6364117) q[3];
sx q[3];
rz(-0.44613871) q[3];
sx q[3];
rz(1.9193468) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196359) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(-1.5420445) q[0];
rz(-2.385335) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(-0.33682987) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23754691) q[0];
sx q[0];
rz(-1.2531452) q[0];
sx q[0];
rz(-1.5730412) q[0];
rz(-pi) q[1];
rz(-2.1977399) q[2];
sx q[2];
rz(-0.34239951) q[2];
sx q[2];
rz(-1.2417718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.017180786) q[1];
sx q[1];
rz(-1.681911) q[1];
sx q[1];
rz(3.0558056) q[1];
x q[2];
rz(0.41542094) q[3];
sx q[3];
rz(-1.0491228) q[3];
sx q[3];
rz(-3.0565302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6280262) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(2.8619518) q[2];
rz(-0.63129342) q[3];
sx q[3];
rz(-2.203233) q[3];
sx q[3];
rz(0.62425557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.30017988) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(-0.77371669) q[1];
sx q[1];
rz(-2.5061889) q[1];
sx q[1];
rz(-2.9255964) q[1];
rz(2.5071267) q[2];
sx q[2];
rz(-0.96426156) q[2];
sx q[2];
rz(-2.5360863) q[2];
rz(1.5935043) q[3];
sx q[3];
rz(-1.1946214) q[3];
sx q[3];
rz(-0.0047636845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
