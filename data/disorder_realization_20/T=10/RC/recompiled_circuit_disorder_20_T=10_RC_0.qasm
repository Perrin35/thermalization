OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.5469172) q[0];
sx q[0];
rz(-2.1407054) q[0];
sx q[0];
rz(0.21240182) q[0];
rz(0.71495932) q[1];
sx q[1];
rz(3.9290805) q[1];
sx q[1];
rz(10.706283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7377388) q[0];
sx q[0];
rz(-1.3699342) q[0];
sx q[0];
rz(1.2234729) q[0];
x q[1];
rz(0.035543156) q[2];
sx q[2];
rz(-2.4487552) q[2];
sx q[2];
rz(2.7685744) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4557588) q[1];
sx q[1];
rz(-1.403192) q[1];
sx q[1];
rz(-1.140825) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1164066) q[3];
sx q[3];
rz(-2.9985399) q[3];
sx q[3];
rz(-1.5798626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.5167351) q[2];
rz(2.8939698) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6973998) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(-2.9852988) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(-1.4888391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.981502) q[0];
sx q[0];
rz(-2.1182051) q[0];
sx q[0];
rz(-1.7574969) q[0];
rz(-pi) q[1];
rz(-2.3180914) q[2];
sx q[2];
rz(-2.245792) q[2];
sx q[2];
rz(1.0456955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9298047) q[1];
sx q[1];
rz(-2.1165407) q[1];
sx q[1];
rz(-2.2499229) q[1];
x q[2];
rz(1.5403455) q[3];
sx q[3];
rz(-1.4342562) q[3];
sx q[3];
rz(-2.2583435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(-2.8448811) q[2];
rz(0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(-1.9077574) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333703) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(2.6112774) q[0];
rz(-0.5439533) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(-1.9525607) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30071242) q[0];
sx q[0];
rz(-1.2749201) q[0];
sx q[0];
rz(1.7353252) q[0];
x q[1];
rz(-1.6109153) q[2];
sx q[2];
rz(-1.3228068) q[2];
sx q[2];
rz(0.2085533) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51057286) q[1];
sx q[1];
rz(-0.90386183) q[1];
sx q[1];
rz(-2.2651947) q[1];
rz(2.6044106) q[3];
sx q[3];
rz(-0.77386412) q[3];
sx q[3];
rz(2.8230132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.8900324) q[2];
rz(-0.45423147) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(-0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(1.9569424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9979447) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(-0.68302897) q[0];
rz(-2.0834288) q[2];
sx q[2];
rz(-1.1626273) q[2];
sx q[2];
rz(0.689091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9222316) q[1];
sx q[1];
rz(-0.9833828) q[1];
sx q[1];
rz(-3.0012793) q[1];
rz(-pi) q[2];
x q[2];
rz(2.437856) q[3];
sx q[3];
rz(-2.6061213) q[3];
sx q[3];
rz(-0.4243917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(0.63956368) q[2];
rz(0.8574287) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(-0.48979428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63067591) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(0.72881126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279014) q[0];
sx q[0];
rz(-1.3662845) q[0];
sx q[0];
rz(-1.6708899) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.38703) q[2];
sx q[2];
rz(-0.73256058) q[2];
sx q[2];
rz(-1.5549321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.049206991) q[1];
sx q[1];
rz(-1.2485463) q[1];
sx q[1];
rz(2.9258123) q[1];
rz(-2.2767931) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(-1.4454696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0430498) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(0.70971242) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(-3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0387886) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(2.2578755) q[0];
rz(-0.9206413) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(-2.9439435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5699428) q[0];
sx q[0];
rz(-0.93402223) q[0];
sx q[0];
rz(-2.3417579) q[0];
x q[1];
rz(-2.501776) q[2];
sx q[2];
rz(-3.0090927) q[2];
sx q[2];
rz(-2.8853531) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3870961) q[1];
sx q[1];
rz(-1.6134141) q[1];
sx q[1];
rz(-2.6864762) q[1];
rz(-pi) q[2];
rz(0.88455172) q[3];
sx q[3];
rz(-0.93730799) q[3];
sx q[3];
rz(1.3214878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(3.0659884) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41912115) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(-1.9627337) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.1694318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5116611) q[0];
sx q[0];
rz(-2.9637664) q[0];
sx q[0];
rz(-1.1849665) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39016907) q[2];
sx q[2];
rz(-1.8346268) q[2];
sx q[2];
rz(-2.2389776) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84050256) q[1];
sx q[1];
rz(-1.2355348) q[1];
sx q[1];
rz(2.5696978) q[1];
x q[2];
rz(-2.7536105) q[3];
sx q[3];
rz(-1.6098607) q[3];
sx q[3];
rz(-0.29400533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.823267) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(0.014483359) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(-1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2248174) q[0];
sx q[0];
rz(-0.48368925) q[0];
sx q[0];
rz(-1.7513357) q[0];
rz(-pi) q[1];
rz(-2.2675603) q[2];
sx q[2];
rz(-1.1296141) q[2];
sx q[2];
rz(2.2639084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22842912) q[1];
sx q[1];
rz(-0.73691165) q[1];
sx q[1];
rz(-1.5494898) q[1];
x q[2];
rz(1.0192972) q[3];
sx q[3];
rz(-0.5987474) q[3];
sx q[3];
rz(2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3056425) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(0.28253728) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95865059) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(-1.4022934) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.8797849) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14784797) q[0];
sx q[0];
rz(-1.9980668) q[0];
sx q[0];
rz(1.8673613) q[0];
rz(-pi) q[1];
rz(-2.8724573) q[2];
sx q[2];
rz(-2.4705774) q[2];
sx q[2];
rz(3.105643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.251235) q[1];
sx q[1];
rz(-2.2489378) q[1];
sx q[1];
rz(-1.8974341) q[1];
rz(-1.3654361) q[3];
sx q[3];
rz(-2.0911651) q[3];
sx q[3];
rz(-0.53232771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.6837439) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(-2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(-1.9816459) q[0];
rz(3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(-1.0704401) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72296952) q[0];
sx q[0];
rz(-1.4040935) q[0];
sx q[0];
rz(-0.82794257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6106597) q[2];
sx q[2];
rz(-0.76518744) q[2];
sx q[2];
rz(2.6305692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66026238) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(0.34923133) q[1];
rz(-pi) q[2];
rz(0.6108547) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(-2.6328997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(1.6147511) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51331818) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(1.0340446) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(-0.037739567) q[2];
sx q[2];
rz(-2.3024617) q[2];
sx q[2];
rz(3.0842177) q[2];
rz(0.74450775) q[3];
sx q[3];
rz(-2.4662938) q[3];
sx q[3];
rz(0.288356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
