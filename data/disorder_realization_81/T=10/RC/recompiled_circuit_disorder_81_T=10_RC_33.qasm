OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(1.3786432) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(-2.690697) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9076951) q[0];
sx q[0];
rz(-1.6595708) q[0];
sx q[0];
rz(1.5560454) q[0];
x q[1];
rz(-2.0499174) q[2];
sx q[2];
rz(-1.4716822) q[2];
sx q[2];
rz(-1.6247768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2921819) q[1];
sx q[1];
rz(-1.0725478) q[1];
sx q[1];
rz(-1.317418) q[1];
rz(1.6879184) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(-2.5205034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(-0.84428865) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-2.5355693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(1.9553604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378368) q[0];
sx q[0];
rz(-0.058996011) q[0];
sx q[0];
rz(-2.8179413) q[0];
rz(1.6329174) q[2];
sx q[2];
rz(-1.7698235) q[2];
sx q[2];
rz(1.3028499) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9895049) q[1];
sx q[1];
rz(-0.67645914) q[1];
sx q[1];
rz(1.3701887) q[1];
x q[2];
rz(0.22963345) q[3];
sx q[3];
rz(-2.4335055) q[3];
sx q[3];
rz(-1.9830444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.9821232) q[2];
rz(0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(2.3348715) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(2.321373) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0778753) q[0];
sx q[0];
rz(-2.2532007) q[0];
sx q[0];
rz(-0.14494411) q[0];
x q[1];
rz(2.9224612) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(0.52106524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8169176) q[1];
sx q[1];
rz(-1.1500689) q[1];
sx q[1];
rz(-2.7540728) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1060171) q[3];
sx q[3];
rz(-1.6764063) q[3];
sx q[3];
rz(-1.6184023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(2.9860949) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(0.43831929) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.320425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5116918) q[0];
sx q[0];
rz(-1.9803847) q[0];
sx q[0];
rz(-3.1303309) q[0];
x q[1];
rz(-2.690372) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(-0.95552432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3418386) q[1];
sx q[1];
rz(-2.2463887) q[1];
sx q[1];
rz(1.9104596) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60453316) q[3];
sx q[3];
rz(-0.88450888) q[3];
sx q[3];
rz(1.3674919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(2.7992115) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8300366) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(-0.86529055) q[0];
rz(-1.226549) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(1.8409761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7428955) q[0];
sx q[0];
rz(-1.5455751) q[0];
sx q[0];
rz(-2.7874649) q[0];
rz(-pi) q[1];
rz(1.498921) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(-0.14771151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3805441) q[1];
sx q[1];
rz(-2.492978) q[1];
sx q[1];
rz(2.5785239) q[1];
rz(-2.8793094) q[3];
sx q[3];
rz(-1.4296921) q[3];
sx q[3];
rz(1.9067681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(0.47719964) q[2];
rz(-2.9495083) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(2.5911962) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5884429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9827305) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(-1.1984675) q[0];
rz(1.6278218) q[2];
sx q[2];
rz(-2.7751623) q[2];
sx q[2];
rz(-2.1232405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99332422) q[1];
sx q[1];
rz(-2.4273708) q[1];
sx q[1];
rz(0.12970129) q[1];
x q[2];
rz(-1.4828959) q[3];
sx q[3];
rz(-2.4368736) q[3];
sx q[3];
rz(0.3375012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(1.7717308) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-1.1184568) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5605374) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1099694) q[0];
sx q[0];
rz(-2.1620746) q[0];
sx q[0];
rz(-1.6615608) q[0];
x q[1];
rz(3.1408429) q[2];
sx q[2];
rz(-3.0095518) q[2];
sx q[2];
rz(-3.0291639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44927412) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(-1.7983789) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3903923) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(-0.62869149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(2.1288669) q[2];
rz(1.9536473) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-2.4429328) q[0];
rz(2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.8922071) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4417393) q[0];
sx q[0];
rz(-0.2158567) q[0];
sx q[0];
rz(2.9237843) q[0];
rz(1.0480568) q[2];
sx q[2];
rz(-1.153423) q[2];
sx q[2];
rz(2.0975031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74832143) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(3.1194411) q[1];
rz(-pi) q[2];
rz(1.1864248) q[3];
sx q[3];
rz(-2.4723408) q[3];
sx q[3];
rz(2.3068908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32968783) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(-1.9630986) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-2.7594574) q[3];
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
rz(1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(-1.4216084) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(-0.59757772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9874728) q[0];
sx q[0];
rz(-0.14053908) q[0];
sx q[0];
rz(1.2921635) q[0];
x q[1];
rz(2.3915646) q[2];
sx q[2];
rz(-0.70493297) q[2];
sx q[2];
rz(-2.7761369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4328879) q[1];
sx q[1];
rz(-1.27379) q[1];
sx q[1];
rz(1.2328641) q[1];
rz(-pi) q[2];
rz(-1.647801) q[3];
sx q[3];
rz(-0.98620755) q[3];
sx q[3];
rz(-1.0494174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9188345) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(2.4694209) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7260872) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(2.0477717) q[0];
rz(-2.6300738) q[2];
sx q[2];
rz(-1.5626972) q[2];
sx q[2];
rz(0.17009232) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2031659) q[1];
sx q[1];
rz(-2.2397579) q[1];
sx q[1];
rz(2.0069564) q[1];
x q[2];
rz(-2.7087595) q[3];
sx q[3];
rz(-1.4308617) q[3];
sx q[3];
rz(2.7452552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(2.771634) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(0.72369408) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-1.205668) q[2];
sx q[2];
rz(-2.7177313) q[2];
sx q[2];
rz(2.360366) q[2];
rz(1.1643812) q[3];
sx q[3];
rz(-1.2953399) q[3];
sx q[3];
rz(-3.0084707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
