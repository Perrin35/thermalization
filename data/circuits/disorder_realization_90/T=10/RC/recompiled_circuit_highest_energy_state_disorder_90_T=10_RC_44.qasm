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
rz(0.69819063) q[0];
sx q[0];
rz(-0.31261045) q[0];
sx q[0];
rz(-0.99553776) q[0];
rz(-2.0350463) q[1];
sx q[1];
rz(-2.3732329) q[1];
sx q[1];
rz(0.33222693) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8790008) q[0];
sx q[0];
rz(-2.342412) q[0];
sx q[0];
rz(-2.5627717) q[0];
rz(-pi) q[1];
rz(1.1071476) q[2];
sx q[2];
rz(-1.6939179) q[2];
sx q[2];
rz(2.5277918) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54094568) q[1];
sx q[1];
rz(-1.9917734) q[1];
sx q[1];
rz(1.0599815) q[1];
rz(-0.99523441) q[3];
sx q[3];
rz(-1.1249295) q[3];
sx q[3];
rz(0.80328548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.535061) q[2];
sx q[2];
rz(-1.0142832) q[2];
sx q[2];
rz(-3.0895341) q[2];
rz(-2.0590797) q[3];
sx q[3];
rz(-2.1942997) q[3];
sx q[3];
rz(1.1516217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.624619) q[0];
sx q[0];
rz(-0.032051429) q[0];
sx q[0];
rz(2.8042941) q[0];
rz(-2.9889122) q[1];
sx q[1];
rz(-2.3744507) q[1];
sx q[1];
rz(-1.5229567) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168287) q[0];
sx q[0];
rz(-1.404002) q[0];
sx q[0];
rz(-2.5186064) q[0];
rz(-0.27045336) q[2];
sx q[2];
rz(-2.4861397) q[2];
sx q[2];
rz(1.356911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26451787) q[1];
sx q[1];
rz(-1.3691069) q[1];
sx q[1];
rz(-1.1796477) q[1];
rz(-pi) q[2];
rz(-2.5378599) q[3];
sx q[3];
rz(-1.4430178) q[3];
sx q[3];
rz(-1.4462245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0308257) q[2];
sx q[2];
rz(-2.775707) q[2];
sx q[2];
rz(-2.7552674) q[2];
rz(-1.857916) q[3];
sx q[3];
rz(-1.176703) q[3];
sx q[3];
rz(-1.2181237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.06685) q[0];
sx q[0];
rz(-1.0800986) q[0];
sx q[0];
rz(2.1225488) q[0];
rz(-0.77541238) q[1];
sx q[1];
rz(-1.8653899) q[1];
sx q[1];
rz(-2.6554241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3306823) q[0];
sx q[0];
rz(-0.031477246) q[0];
sx q[0];
rz(-3.0481008) q[0];
x q[1];
rz(2.0569226) q[2];
sx q[2];
rz(-2.297431) q[2];
sx q[2];
rz(-1.655533) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61709014) q[1];
sx q[1];
rz(-0.5382584) q[1];
sx q[1];
rz(-2.7048802) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8649549) q[3];
sx q[3];
rz(-1.8801062) q[3];
sx q[3];
rz(-0.33634392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7458618) q[2];
sx q[2];
rz(-1.9755325) q[2];
sx q[2];
rz(-2.3365848) q[2];
rz(-2.4896367) q[3];
sx q[3];
rz(-1.1019573) q[3];
sx q[3];
rz(0.93418795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.4790799) q[0];
sx q[0];
rz(-1.3211687) q[0];
sx q[0];
rz(1.3385734) q[0];
rz(1.592912) q[1];
sx q[1];
rz(-0.70392307) q[1];
sx q[1];
rz(2.7159363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073624728) q[0];
sx q[0];
rz(-1.5626161) q[0];
sx q[0];
rz(0.091853022) q[0];
x q[1];
rz(1.2629444) q[2];
sx q[2];
rz(-1.6037707) q[2];
sx q[2];
rz(0.82885447) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7158) q[1];
sx q[1];
rz(-1.1669901) q[1];
sx q[1];
rz(0.76622643) q[1];
rz(-pi) q[2];
x q[2];
rz(1.852774) q[3];
sx q[3];
rz(-2.2681142) q[3];
sx q[3];
rz(-0.048843637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8694596) q[2];
sx q[2];
rz(-1.3674066) q[2];
sx q[2];
rz(-2.7222705) q[2];
rz(1.9648633) q[3];
sx q[3];
rz(-0.71507016) q[3];
sx q[3];
rz(-0.56593219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6856573) q[0];
sx q[0];
rz(-2.3929907) q[0];
sx q[0];
rz(-1.5022044) q[0];
rz(-1.7078851) q[1];
sx q[1];
rz(-1.8610443) q[1];
sx q[1];
rz(2.7388403) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9929189) q[0];
sx q[0];
rz(-0.67280992) q[0];
sx q[0];
rz(1.1514173) q[0];
x q[1];
rz(-1.5115159) q[2];
sx q[2];
rz(-1.0067954) q[2];
sx q[2];
rz(0.21280542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50583831) q[1];
sx q[1];
rz(-0.16316667) q[1];
sx q[1];
rz(-2.5058305) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4536269) q[3];
sx q[3];
rz(-2.2343383) q[3];
sx q[3];
rz(2.508916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5002284) q[2];
sx q[2];
rz(-2.0945175) q[2];
sx q[2];
rz(0.25263986) q[2];
rz(-2.3371475) q[3];
sx q[3];
rz(-0.52144709) q[3];
sx q[3];
rz(-0.49853244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8135524) q[0];
sx q[0];
rz(-0.10843065) q[0];
sx q[0];
rz(-2.257708) q[0];
rz(0.49993316) q[1];
sx q[1];
rz(-1.6729313) q[1];
sx q[1];
rz(0.51441851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2470236) q[0];
sx q[0];
rz(-1.2366364) q[0];
sx q[0];
rz(0.96118013) q[0];
rz(1.8158004) q[2];
sx q[2];
rz(-1.0520237) q[2];
sx q[2];
rz(3.0381448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0658256) q[1];
sx q[1];
rz(-1.3525241) q[1];
sx q[1];
rz(1.1595919) q[1];
x q[2];
rz(-1.5443956) q[3];
sx q[3];
rz(-2.3599778) q[3];
sx q[3];
rz(-1.7780245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9875235) q[2];
sx q[2];
rz(-1.6776626) q[2];
sx q[2];
rz(-0.12518159) q[2];
rz(2.3212738) q[3];
sx q[3];
rz(-1.9671974) q[3];
sx q[3];
rz(-2.7952747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5888551) q[0];
sx q[0];
rz(-0.6468361) q[0];
sx q[0];
rz(3.0297025) q[0];
rz(-2.0018068) q[1];
sx q[1];
rz(-0.21509376) q[1];
sx q[1];
rz(1.8591759) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70045602) q[0];
sx q[0];
rz(-2.4917291) q[0];
sx q[0];
rz(2.031424) q[0];
rz(-pi) q[1];
rz(-1.4256023) q[2];
sx q[2];
rz(-1.34158) q[2];
sx q[2];
rz(-0.21266937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0874028) q[1];
sx q[1];
rz(-0.8381745) q[1];
sx q[1];
rz(-2.0778632) q[1];
rz(-pi) q[2];
rz(-2.4573648) q[3];
sx q[3];
rz(-1.5260395) q[3];
sx q[3];
rz(2.0852722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.68975) q[2];
sx q[2];
rz(-0.30561438) q[2];
sx q[2];
rz(-1.8243194) q[2];
rz(-0.02056038) q[3];
sx q[3];
rz(-2.8097184) q[3];
sx q[3];
rz(2.1338972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2059712) q[0];
sx q[0];
rz(-0.99140778) q[0];
sx q[0];
rz(0.29888612) q[0];
rz(2.0805953) q[1];
sx q[1];
rz(-1.3875049) q[1];
sx q[1];
rz(-1.908173) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3372634) q[0];
sx q[0];
rz(-1.6676039) q[0];
sx q[0];
rz(1.2362572) q[0];
rz(1.226108) q[2];
sx q[2];
rz(-0.12154254) q[2];
sx q[2];
rz(-0.28608382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2005991) q[1];
sx q[1];
rz(-0.72807136) q[1];
sx q[1];
rz(2.0242101) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0954433) q[3];
sx q[3];
rz(-0.97844175) q[3];
sx q[3];
rz(0.86529675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79248205) q[2];
sx q[2];
rz(-1.0976378) q[2];
sx q[2];
rz(2.9343572) q[2];
rz(-0.66050291) q[3];
sx q[3];
rz(-2.1410172) q[3];
sx q[3];
rz(-1.8628666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58885634) q[0];
sx q[0];
rz(-1.4016466) q[0];
sx q[0];
rz(-1.9197585) q[0];
rz(2.6264722) q[1];
sx q[1];
rz(-3.0312067) q[1];
sx q[1];
rz(2.5882904) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.57896) q[0];
sx q[0];
rz(-1.4427358) q[0];
sx q[0];
rz(2.1306031) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9917351) q[2];
sx q[2];
rz(-1.6302262) q[2];
sx q[2];
rz(0.42039117) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8749145) q[1];
sx q[1];
rz(-1.5767225) q[1];
sx q[1];
rz(1.3937217) q[1];
rz(-pi) q[2];
rz(-2.4724602) q[3];
sx q[3];
rz(-1.5894984) q[3];
sx q[3];
rz(1.6220055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33636609) q[2];
sx q[2];
rz(-1.1104501) q[2];
sx q[2];
rz(2.6160348) q[2];
rz(-1.4793652) q[3];
sx q[3];
rz(-0.95366228) q[3];
sx q[3];
rz(-2.3814538) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084759921) q[0];
sx q[0];
rz(-0.79462093) q[0];
sx q[0];
rz(0.42924616) q[0];
rz(-1.5224573) q[1];
sx q[1];
rz(-1.2702962) q[1];
sx q[1];
rz(-0.92990184) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36741716) q[0];
sx q[0];
rz(-0.40433592) q[0];
sx q[0];
rz(1.0003759) q[0];
rz(-1.11082) q[2];
sx q[2];
rz(-0.24402741) q[2];
sx q[2];
rz(-2.8619253) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4082171) q[1];
sx q[1];
rz(-1.2885417) q[1];
sx q[1];
rz(0.1677081) q[1];
rz(1.6305805) q[3];
sx q[3];
rz(-1.0185223) q[3];
sx q[3];
rz(-2.5552101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.736019) q[2];
sx q[2];
rz(-2.6655727) q[2];
sx q[2];
rz(-2.948577) q[2];
rz(0.080282601) q[3];
sx q[3];
rz(-2.2716227) q[3];
sx q[3];
rz(-1.75753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3748462) q[0];
sx q[0];
rz(-1.9369047) q[0];
sx q[0];
rz(1.9932224) q[0];
rz(-2.171352) q[1];
sx q[1];
rz(-1.932825) q[1];
sx q[1];
rz(1.8995151) q[1];
rz(-1.8850897) q[2];
sx q[2];
rz(-2.3705924) q[2];
sx q[2];
rz(2.3938897) q[2];
rz(1.5469848) q[3];
sx q[3];
rz(-1.4766558) q[3];
sx q[3];
rz(-2.7500931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
