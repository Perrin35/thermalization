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
rz(-0.22696683) q[0];
sx q[0];
rz(-2.5300955) q[0];
sx q[0];
rz(-0.55846941) q[0];
rz(-2.5418169) q[1];
sx q[1];
rz(-1.37473) q[1];
sx q[1];
rz(3.0629646) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0538875) q[0];
sx q[0];
rz(-2.1057406) q[0];
sx q[0];
rz(-0.91607262) q[0];
x q[1];
rz(-2.4628381) q[2];
sx q[2];
rz(-1.6683443) q[2];
sx q[2];
rz(0.3435979) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1153733) q[1];
sx q[1];
rz(-1.3182606) q[1];
sx q[1];
rz(1.4863243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8923678) q[3];
sx q[3];
rz(-1.235629) q[3];
sx q[3];
rz(-0.17091076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2416396) q[2];
sx q[2];
rz(-3.0934379) q[2];
sx q[2];
rz(-1.0345577) q[2];
rz(3.0311846) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(-0.30606562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0050874) q[0];
sx q[0];
rz(-1.6965447) q[0];
sx q[0];
rz(1.9553631) q[0];
rz(-1.2677445) q[1];
sx q[1];
rz(-0.45919752) q[1];
sx q[1];
rz(-1.3892106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25042501) q[0];
sx q[0];
rz(-1.0396663) q[0];
sx q[0];
rz(2.7851339) q[0];
rz(-pi) q[1];
rz(2.9302338) q[2];
sx q[2];
rz(-0.49079259) q[2];
sx q[2];
rz(-0.15891128) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0276549) q[1];
sx q[1];
rz(-1.9009556) q[1];
sx q[1];
rz(0.3056674) q[1];
rz(0.94654437) q[3];
sx q[3];
rz(-1.6485751) q[3];
sx q[3];
rz(-2.3619224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9691951) q[2];
sx q[2];
rz(-2.0266504) q[2];
sx q[2];
rz(-1.4364852) q[2];
rz(-1.8819594) q[3];
sx q[3];
rz(-2.6862222) q[3];
sx q[3];
rz(-0.026738515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7773975) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(-1.9245603) q[0];
rz(-2.9265535) q[1];
sx q[1];
rz(-1.2478849) q[1];
sx q[1];
rz(-1.4395813) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3830945) q[0];
sx q[0];
rz(-1.4466084) q[0];
sx q[0];
rz(0.66207768) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0820175) q[2];
sx q[2];
rz(-1.4892092) q[2];
sx q[2];
rz(-2.4961584) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7172723) q[1];
sx q[1];
rz(-1.3417146) q[1];
sx q[1];
rz(0.13026625) q[1];
rz(-pi) q[2];
rz(-1.2974627) q[3];
sx q[3];
rz(-1.7810827) q[3];
sx q[3];
rz(-3.1137636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3415459) q[2];
sx q[2];
rz(-1.6365106) q[2];
sx q[2];
rz(2.3254507) q[2];
rz(0.0084361313) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(1.5754023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4207882) q[0];
sx q[0];
rz(-1.9705462) q[0];
sx q[0];
rz(-2.8724443) q[0];
rz(1.7587657) q[1];
sx q[1];
rz(-0.56126422) q[1];
sx q[1];
rz(-0.47163481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1203493) q[0];
sx q[0];
rz(-0.85130748) q[0];
sx q[0];
rz(-2.9553652) q[0];
x q[1];
rz(3.0582171) q[2];
sx q[2];
rz(-0.40503392) q[2];
sx q[2];
rz(1.1277367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9593352) q[1];
sx q[1];
rz(-1.8885837) q[1];
sx q[1];
rz(2.5885568) q[1];
x q[2];
rz(-2.7876623) q[3];
sx q[3];
rz(-2.2056863) q[3];
sx q[3];
rz(-0.77087444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4857594) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(-1.887623) q[2];
rz(-2.8433825) q[3];
sx q[3];
rz(-1.2757653) q[3];
sx q[3];
rz(1.5378753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0373847) q[0];
sx q[0];
rz(-0.91008121) q[0];
sx q[0];
rz(0.75697672) q[0];
rz(1.0589927) q[1];
sx q[1];
rz(-1.6964361) q[1];
sx q[1];
rz(-0.34245488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7624906) q[0];
sx q[0];
rz(-2.6978793) q[0];
sx q[0];
rz(2.0279473) q[0];
rz(-pi) q[1];
rz(-0.86910291) q[2];
sx q[2];
rz(-1.822374) q[2];
sx q[2];
rz(-0.053843018) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96794242) q[1];
sx q[1];
rz(-1.0999803) q[1];
sx q[1];
rz(-1.985926) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60635546) q[3];
sx q[3];
rz(-2.6487051) q[3];
sx q[3];
rz(2.3541401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4464438) q[2];
sx q[2];
rz(-1.4651639) q[2];
sx q[2];
rz(-0.62015074) q[2];
rz(0.54689637) q[3];
sx q[3];
rz(-0.20709012) q[3];
sx q[3];
rz(2.0641522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3243489) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(1.9203583) q[0];
rz(-2.0381894) q[1];
sx q[1];
rz(-1.9788479) q[1];
sx q[1];
rz(0.81072909) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4909523) q[0];
sx q[0];
rz(-1.8260667) q[0];
sx q[0];
rz(1.4064691) q[0];
rz(-pi) q[1];
rz(2.4898363) q[2];
sx q[2];
rz(-2.0988191) q[2];
sx q[2];
rz(-0.17490696) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.74720464) q[1];
sx q[1];
rz(-1.516368) q[1];
sx q[1];
rz(-1.4318951) q[1];
rz(2.8437988) q[3];
sx q[3];
rz(-1.9341365) q[3];
sx q[3];
rz(0.17223528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3409884) q[2];
sx q[2];
rz(-1.1977414) q[2];
sx q[2];
rz(1.9761337) q[2];
rz(-3.0321339) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(0.060733184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8282181) q[0];
sx q[0];
rz(-0.010638588) q[0];
sx q[0];
rz(0.64205545) q[0];
rz(-2.8984046) q[1];
sx q[1];
rz(-0.41047341) q[1];
sx q[1];
rz(-0.64116716) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12578861) q[0];
sx q[0];
rz(-2.1282853) q[0];
sx q[0];
rz(2.1899738) q[0];
rz(1.2095318) q[2];
sx q[2];
rz(-2.2881977) q[2];
sx q[2];
rz(-0.016005767) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4879681) q[1];
sx q[1];
rz(-1.2418126) q[1];
sx q[1];
rz(2.3423561) q[1];
rz(2.9092714) q[3];
sx q[3];
rz(-1.9594176) q[3];
sx q[3];
rz(-0.53750932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0082561) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(-2.1978417) q[2];
rz(0.64822316) q[3];
sx q[3];
rz(-2.4002878) q[3];
sx q[3];
rz(0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7159395) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(-2.7124523) q[0];
rz(-0.40402135) q[1];
sx q[1];
rz(-0.41530135) q[1];
sx q[1];
rz(-0.9187575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2935393) q[0];
sx q[0];
rz(-1.812879) q[0];
sx q[0];
rz(-2.5771067) q[0];
x q[1];
rz(1.3098251) q[2];
sx q[2];
rz(-1.6918285) q[2];
sx q[2];
rz(-2.3132035) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1240108) q[1];
sx q[1];
rz(-0.38330829) q[1];
sx q[1];
rz(-0.95335754) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2567344) q[3];
sx q[3];
rz(-1.6410284) q[3];
sx q[3];
rz(-2.8633871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.018844133) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(1.7913691) q[2];
rz(2.8090779) q[3];
sx q[3];
rz(-0.88034383) q[3];
sx q[3];
rz(-1.3688709) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3933082) q[0];
sx q[0];
rz(-2.4916593) q[0];
sx q[0];
rz(-2.6848324) q[0];
rz(-1.4211753) q[1];
sx q[1];
rz(-1.8616118) q[1];
sx q[1];
rz(-3.0866887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0676014) q[0];
sx q[0];
rz(-1.5996063) q[0];
sx q[0];
rz(-2.3037698) q[0];
x q[1];
rz(-2.0012399) q[2];
sx q[2];
rz(-1.9483742) q[2];
sx q[2];
rz(2.0098639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0662288) q[1];
sx q[1];
rz(-1.3540596) q[1];
sx q[1];
rz(-1.073887) q[1];
rz(-pi) q[2];
rz(-2.6056882) q[3];
sx q[3];
rz(-1.2782989) q[3];
sx q[3];
rz(1.747594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88825893) q[2];
sx q[2];
rz(-0.71037018) q[2];
sx q[2];
rz(2.9448275) q[2];
rz(-2.0678068) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(2.4055068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1573023) q[0];
sx q[0];
rz(-0.8534011) q[0];
sx q[0];
rz(-0.89944696) q[0];
rz(2.8304214) q[1];
sx q[1];
rz(-1.9322461) q[1];
sx q[1];
rz(1.4003632) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3547826) q[0];
sx q[0];
rz(-1.3616832) q[0];
sx q[0];
rz(-2.1378921) q[0];
rz(-1.1323523) q[2];
sx q[2];
rz(-0.52784) q[2];
sx q[2];
rz(2.6486625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8099762) q[1];
sx q[1];
rz(-0.73996937) q[1];
sx q[1];
rz(-2.5472067) q[1];
rz(-pi) q[2];
rz(0.71615852) q[3];
sx q[3];
rz(-2.0422438) q[3];
sx q[3];
rz(0.52385073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.060999) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(-1.9798123) q[2];
rz(-1.0787841) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(0.88206464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5494736) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(0.1334162) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(1.5897122) q[2];
sx q[2];
rz(-3.0881791) q[2];
sx q[2];
rz(-1.5759158) q[2];
rz(0.50167636) q[3];
sx q[3];
rz(-1.5080843) q[3];
sx q[3];
rz(2.9894473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
