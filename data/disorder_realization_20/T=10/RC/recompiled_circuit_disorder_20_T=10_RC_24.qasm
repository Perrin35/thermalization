OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(-0.21240182) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(1.8600872) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6706657) q[0];
sx q[0];
rz(-2.7424194) q[0];
sx q[0];
rz(-1.0317208) q[0];
rz(-pi) q[1];
rz(-0.035543156) q[2];
sx q[2];
rz(-2.4487552) q[2];
sx q[2];
rz(-2.7685744) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3329831) q[1];
sx q[1];
rz(-1.994351) q[1];
sx q[1];
rz(2.9575648) q[1];
rz(-pi) q[2];
rz(3.078457) q[3];
sx q[3];
rz(-1.4423443) q[3];
sx q[3];
rz(-2.0383143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(1.6248576) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(2.9852988) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(1.6527536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50870897) q[0];
sx q[0];
rz(-1.4116305) q[0];
sx q[0];
rz(0.55523086) q[0];
x q[1];
rz(2.4376051) q[2];
sx q[2];
rz(-2.1805602) q[2];
sx q[2];
rz(-0.068254452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3540346) q[1];
sx q[1];
rz(-2.2984142) q[1];
sx q[1];
rz(-2.3393199) q[1];
rz(3.00499) q[3];
sx q[3];
rz(-1.540629) q[3];
sx q[3];
rz(-2.4581916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8504101) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(-2.8448811) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(-1.2338352) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082224) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(-2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.9525607) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2169164) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-0.49305537) q[0];
x q[1];
rz(0.24818111) q[2];
sx q[2];
rz(-1.5319053) q[2];
sx q[2];
rz(1.789202) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4424858) q[1];
sx q[1];
rz(-2.2190296) q[1];
sx q[1];
rz(-0.68251619) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6044106) q[3];
sx q[3];
rz(-0.77386412) q[3];
sx q[3];
rz(-0.31857946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.2515602) q[2];
rz(-0.45423147) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(-2.5355133) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(-1.1846503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957164) q[0];
sx q[0];
rz(-1.3564975) q[0];
sx q[0];
rz(1.3927668) q[0];
rz(-pi) q[1];
rz(-2.6809533) q[2];
sx q[2];
rz(-2.0377635) q[2];
sx q[2];
rz(0.66191445) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6726917) q[1];
sx q[1];
rz(-2.5395782) q[1];
sx q[1];
rz(-1.7778346) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7167927) q[3];
sx q[3];
rz(-1.9072755) q[3];
sx q[3];
rz(2.6257023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23345315) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(2.502029) q[2];
rz(2.284164) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-0.48979428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(-2.6211991) q[0];
rz(3.0386472) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(0.72881126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8048808) q[0];
sx q[0];
rz(-1.4727955) q[0];
sx q[0];
rz(2.9360807) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9786733) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(1.7999072) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.049206991) q[1];
sx q[1];
rz(-1.2485463) q[1];
sx q[1];
rz(0.21578034) q[1];
rz(-2.2767931) q[3];
sx q[3];
rz(-0.74665657) q[3];
sx q[3];
rz(1.4454696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(2.4318802) q[2];
rz(-2.5085311) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(0.13171296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028041) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(-2.2578755) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(-2.9439435) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55006856) q[0];
sx q[0];
rz(-0.95614377) q[0];
sx q[0];
rz(-0.7556677) q[0];
rz(-pi) q[1];
rz(2.501776) q[2];
sx q[2];
rz(-3.0090927) q[2];
sx q[2];
rz(-0.25623955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27053988) q[1];
sx q[1];
rz(-2.6846243) q[1];
sx q[1];
rz(0.096710042) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.382155) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(0.2021377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(-0.075604288) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(-0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41912115) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(-1.9627337) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(-1.9721608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1202576) q[0];
sx q[0];
rz(-1.7354256) q[0];
sx q[0];
rz(-0.067532587) q[0];
rz(0.39016907) q[2];
sx q[2];
rz(-1.8346268) q[2];
sx q[2];
rz(2.2389776) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6199854) q[1];
sx q[1];
rz(-1.0343401) q[1];
sx q[1];
rz(-1.9636088) q[1];
rz(-pi) q[2];
rz(-0.38798214) q[3];
sx q[3];
rz(-1.6098607) q[3];
sx q[3];
rz(0.29400533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.823267) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(3.1271093) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-0.55317944) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(-1.4315804) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(2.6838578) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2248174) q[0];
sx q[0];
rz(-2.6579034) q[0];
sx q[0];
rz(-1.3902569) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2051815) q[2];
sx q[2];
rz(-0.80447703) q[2];
sx q[2];
rz(-0.22125439) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.22842912) q[1];
sx q[1];
rz(-0.73691165) q[1];
sx q[1];
rz(-1.5921028) q[1];
rz(-pi) q[2];
rz(-2.1222955) q[3];
sx q[3];
rz(-0.5987474) q[3];
sx q[3];
rz(2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3056425) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(2.8590554) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95865059) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.8797849) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5926873) q[0];
sx q[0];
rz(-1.840001) q[0];
sx q[0];
rz(2.6972428) q[0];
rz(-pi) q[1];
rz(-2.8724573) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(-3.105643) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.251235) q[1];
sx q[1];
rz(-2.2489378) q[1];
sx q[1];
rz(-1.8974341) q[1];
rz(-pi) q[2];
rz(0.52957876) q[3];
sx q[3];
rz(-1.3929318) q[3];
sx q[3];
rz(1.1416658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.4578488) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.89501971) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(2.0711526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150072) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(1.8146145) q[0];
x q[1];
rz(-0.038254914) q[2];
sx q[2];
rz(-0.80637156) q[2];
sx q[2];
rz(2.5753266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66214661) q[1];
sx q[1];
rz(-1.3227191) q[1];
sx q[1];
rz(2.3723888) q[1];
x q[2];
rz(1.3954193) q[3];
sx q[3];
rz(-1.8150107) q[3];
sx q[3];
rz(-2.0004686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76876172) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(-1.5268415) q[2];
rz(-0.57957831) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6282745) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(-1.0340446) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(0.8387768) q[2];
sx q[2];
rz(-1.5427187) q[2];
sx q[2];
rz(-1.6533921) q[2];
rz(2.6092929) q[3];
sx q[3];
rz(-1.1333864) q[3];
sx q[3];
rz(-0.65896853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
