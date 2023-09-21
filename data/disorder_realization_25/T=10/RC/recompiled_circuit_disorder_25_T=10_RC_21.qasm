OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(-0.087892428) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9633006) q[0];
sx q[0];
rz(-2.0110197) q[0];
sx q[0];
rz(-2.630169) q[0];
rz(-0.90859969) q[2];
sx q[2];
rz(-2.5052921) q[2];
sx q[2];
rz(-1.5696021) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80884113) q[1];
sx q[1];
rz(-2.07507) q[1];
sx q[1];
rz(1.5732952) q[1];
x q[2];
rz(-1.5486693) q[3];
sx q[3];
rz(-0.25875729) q[3];
sx q[3];
rz(0.72507897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(1.6128444) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(-1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5434718) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.233261) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20576661) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(0.68769023) q[0];
x q[1];
rz(2.9821175) q[2];
sx q[2];
rz(-0.6808241) q[2];
sx q[2];
rz(-2.8603539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.711449) q[1];
sx q[1];
rz(-1.5866382) q[1];
sx q[1];
rz(1.7769017) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66744653) q[3];
sx q[3];
rz(-1.3818041) q[3];
sx q[3];
rz(-2.8611956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1374986) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(0.0022350524) q[2];
rz(2.3114752) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(-1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925791) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-2.2900443) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(3.070389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035633798) q[0];
sx q[0];
rz(-1.7674315) q[0];
sx q[0];
rz(-1.8133624) q[0];
rz(-1.115828) q[2];
sx q[2];
rz(-1.3340545) q[2];
sx q[2];
rz(-2.4617705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60656602) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(-1.4024629) q[1];
rz(2.381388) q[3];
sx q[3];
rz(-1.01902) q[3];
sx q[3];
rz(2.9387568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(-0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-2.6383242) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4098542) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(1.1666974) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(-1.694214) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7839514) q[0];
sx q[0];
rz(-0.28342993) q[0];
sx q[0];
rz(-0.95143239) q[0];
rz(-pi) q[1];
rz(3.1324603) q[2];
sx q[2];
rz(-1.7491241) q[2];
sx q[2];
rz(2.669131) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4059658) q[1];
sx q[1];
rz(-1.4117129) q[1];
sx q[1];
rz(0.87079485) q[1];
rz(2.5656384) q[3];
sx q[3];
rz(-2.7154657) q[3];
sx q[3];
rz(1.061561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(2.9966667) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1458364) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(-2.4575535) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(-1.9285944) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1319259) q[0];
sx q[0];
rz(-1.3103232) q[0];
sx q[0];
rz(-2.7564604) q[0];
rz(0.80981363) q[2];
sx q[2];
rz(-0.81293101) q[2];
sx q[2];
rz(-0.44250689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4124406) q[1];
sx q[1];
rz(-1.8262595) q[1];
sx q[1];
rz(-1.6925807) q[1];
rz(2.9641987) q[3];
sx q[3];
rz(-1.3131485) q[3];
sx q[3];
rz(-2.859476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6325536) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-2.5693494) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260139) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(2.8552326) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(0.58475959) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6385348) q[0];
sx q[0];
rz(-0.98031822) q[0];
sx q[0];
rz(1.1919341) q[0];
rz(-pi) q[1];
rz(3.1044664) q[2];
sx q[2];
rz(-0.6147487) q[2];
sx q[2];
rz(2.6908713) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4510182) q[1];
sx q[1];
rz(-1.8585397) q[1];
sx q[1];
rz(-0.87177999) q[1];
rz(-2.8486409) q[3];
sx q[3];
rz(-2.9679899) q[3];
sx q[3];
rz(-2.0041182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16453234) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(-2.1785054) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(2.5928296) q[0];
rz(3.0896297) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(2.5351977) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604027) q[0];
sx q[0];
rz(-0.58433825) q[0];
sx q[0];
rz(-2.0085745) q[0];
rz(-pi) q[1];
rz(0.90784448) q[2];
sx q[2];
rz(-1.8722476) q[2];
sx q[2];
rz(2.9482258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7402621) q[1];
sx q[1];
rz(-1.4583424) q[1];
sx q[1];
rz(-2.3349891) q[1];
x q[2];
rz(-1.9803195) q[3];
sx q[3];
rz(-1.8990574) q[3];
sx q[3];
rz(1.8144153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-2.6331804) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(-1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-0.34582129) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-0.22769134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84275093) q[0];
sx q[0];
rz(-1.1584873) q[0];
sx q[0];
rz(0.077667872) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.20485) q[2];
sx q[2];
rz(-1.1507251) q[2];
sx q[2];
rz(0.36047381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7714872) q[1];
sx q[1];
rz(-0.58413726) q[1];
sx q[1];
rz(0.41569709) q[1];
rz(3.0759096) q[3];
sx q[3];
rz(-0.55995299) q[3];
sx q[3];
rz(2.5064859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(-0.36455425) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39695981) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(0.80378419) q[0];
rz(-2.0869758) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(-1.1484336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6156857) q[0];
sx q[0];
rz(-1.6427759) q[0];
sx q[0];
rz(3.1176438) q[0];
rz(-0.36275136) q[2];
sx q[2];
rz(-1.0748378) q[2];
sx q[2];
rz(-0.19778684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.25363898) q[1];
sx q[1];
rz(-1.9983074) q[1];
sx q[1];
rz(2.0931787) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6059444) q[3];
sx q[3];
rz(-1.9225549) q[3];
sx q[3];
rz(2.7621321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-1.0174948) q[2];
sx q[2];
rz(0.77825528) q[2];
rz(0.19566472) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739928) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(0.69828066) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2525576) q[0];
sx q[0];
rz(-1.4989984) q[0];
sx q[0];
rz(-0.30307146) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.024784879) q[2];
sx q[2];
rz(-1.6449271) q[2];
sx q[2];
rz(-2.5075846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5729534) q[1];
sx q[1];
rz(-2.0268974) q[1];
sx q[1];
rz(2.9005463) q[1];
rz(-pi) q[2];
rz(-1.2020338) q[3];
sx q[3];
rz(-0.59951111) q[3];
sx q[3];
rz(-2.5332019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0239821) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(-2.3013766) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8828076) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-1.9188948) q[2];
sx q[2];
rz(-0.99708996) q[2];
sx q[2];
rz(-0.70380824) q[2];
rz(2.9914231) q[3];
sx q[3];
rz(-0.88610813) q[3];
sx q[3];
rz(-0.025972493) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];