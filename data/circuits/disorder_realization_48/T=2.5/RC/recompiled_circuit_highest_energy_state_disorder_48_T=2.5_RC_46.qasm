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
rz(2.7895522) q[0];
sx q[0];
rz(-2.3166603) q[0];
sx q[0];
rz(2.6068249) q[0];
rz(-2.1575902) q[1];
sx q[1];
rz(-0.50552955) q[1];
sx q[1];
rz(-1.2619789) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5231294) q[0];
sx q[0];
rz(-0.97181706) q[0];
sx q[0];
rz(1.8904786) q[0];
rz(-pi) q[1];
rz(0.87902576) q[2];
sx q[2];
rz(-2.3874385) q[2];
sx q[2];
rz(1.1471495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.64691862) q[1];
sx q[1];
rz(-0.54975408) q[1];
sx q[1];
rz(0.099262909) q[1];
rz(2.7522911) q[3];
sx q[3];
rz(-1.7598146) q[3];
sx q[3];
rz(2.1019328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24903211) q[2];
sx q[2];
rz(-2.6292215) q[2];
sx q[2];
rz(-1.4830291) q[2];
rz(0.28462166) q[3];
sx q[3];
rz(-2.5183545) q[3];
sx q[3];
rz(2.0774138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677218) q[0];
sx q[0];
rz(-0.72672788) q[0];
sx q[0];
rz(-3.100585) q[0];
rz(-1.9339336) q[1];
sx q[1];
rz(-0.227808) q[1];
sx q[1];
rz(1.4606732) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631993) q[0];
sx q[0];
rz(-1.5413056) q[0];
sx q[0];
rz(0.14568744) q[0];
rz(-2.1182437) q[2];
sx q[2];
rz(-2.259575) q[2];
sx q[2];
rz(2.7163986) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41439287) q[1];
sx q[1];
rz(-1.314394) q[1];
sx q[1];
rz(-3.1016769) q[1];
rz(-1.8051265) q[3];
sx q[3];
rz(-0.70584345) q[3];
sx q[3];
rz(-2.2587551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66551208) q[2];
sx q[2];
rz(-1.6802695) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(-3.0604073) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(1.4359052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.3944893) q[0];
sx q[0];
rz(-3.1293226) q[0];
sx q[0];
rz(-0.61007208) q[0];
rz(-1.3129781) q[1];
sx q[1];
rz(-2.4459631) q[1];
sx q[1];
rz(2.5909766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.074215502) q[0];
sx q[0];
rz(-2.2226637) q[0];
sx q[0];
rz(-1.1135191) q[0];
rz(2.8189893) q[2];
sx q[2];
rz(-1.3122953) q[2];
sx q[2];
rz(0.58835627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50013257) q[1];
sx q[1];
rz(-1.8696897) q[1];
sx q[1];
rz(2.9445093) q[1];
rz(-pi) q[2];
rz(2.9987644) q[3];
sx q[3];
rz(-0.83072829) q[3];
sx q[3];
rz(-2.9252441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6637471) q[2];
sx q[2];
rz(-0.17243324) q[2];
sx q[2];
rz(2.0752068) q[2];
rz(-1.1586698) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(0.042044736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.31024194) q[0];
sx q[0];
rz(-2.5284335) q[0];
sx q[0];
rz(0.9170652) q[0];
rz(0.65545583) q[1];
sx q[1];
rz(-2.2323445) q[1];
sx q[1];
rz(2.7520032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87990078) q[0];
sx q[0];
rz(-1.3845056) q[0];
sx q[0];
rz(-3.0584719) q[0];
rz(-pi) q[1];
rz(-2.5823103) q[2];
sx q[2];
rz(-0.96508316) q[2];
sx q[2];
rz(-2.9756211) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.37871088) q[1];
sx q[1];
rz(-0.40232752) q[1];
sx q[1];
rz(1.4961924) q[1];
rz(-pi) q[2];
rz(0.54549952) q[3];
sx q[3];
rz(-1.0441458) q[3];
sx q[3];
rz(-2.5712476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3098844) q[2];
sx q[2];
rz(-0.94253057) q[2];
sx q[2];
rz(-1.3523098) q[2];
rz(0.74812198) q[3];
sx q[3];
rz(-1.3040521) q[3];
sx q[3];
rz(-1.5266533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8378976) q[0];
sx q[0];
rz(-0.56981531) q[0];
sx q[0];
rz(1.3414398) q[0];
rz(-0.98126423) q[1];
sx q[1];
rz(-1.860268) q[1];
sx q[1];
rz(-0.70995465) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2316596) q[0];
sx q[0];
rz(-1.2450448) q[0];
sx q[0];
rz(-1.5909821) q[0];
x q[1];
rz(2.9030062) q[2];
sx q[2];
rz(-2.1338531) q[2];
sx q[2];
rz(-0.67367879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73588307) q[1];
sx q[1];
rz(-2.8676979) q[1];
sx q[1];
rz(-0.89400684) q[1];
rz(1.4758797) q[3];
sx q[3];
rz(-2.1664985) q[3];
sx q[3];
rz(1.9323521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95825163) q[2];
sx q[2];
rz(-0.7603344) q[2];
sx q[2];
rz(-0.90671268) q[2];
rz(0.18236154) q[3];
sx q[3];
rz(-1.4401108) q[3];
sx q[3];
rz(2.2831634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3458503) q[0];
sx q[0];
rz(-0.95473552) q[0];
sx q[0];
rz(-2.9398651) q[0];
rz(1.7851104) q[1];
sx q[1];
rz(-0.67142612) q[1];
sx q[1];
rz(1.1511525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042554458) q[0];
sx q[0];
rz(-1.7793613) q[0];
sx q[0];
rz(-2.0046356) q[0];
rz(2.5348659) q[2];
sx q[2];
rz(-1.1836241) q[2];
sx q[2];
rz(-0.78578709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6342625) q[1];
sx q[1];
rz(-1.9729904) q[1];
sx q[1];
rz(-0.62437727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3428976) q[3];
sx q[3];
rz(-1.0712475) q[3];
sx q[3];
rz(-1.4394906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98728937) q[2];
sx q[2];
rz(-2.5550948) q[2];
sx q[2];
rz(-2.7941373) q[2];
rz(-0.82184982) q[3];
sx q[3];
rz(-1.776266) q[3];
sx q[3];
rz(-1.52389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3407985) q[0];
sx q[0];
rz(-2.4947385) q[0];
sx q[0];
rz(2.0183753) q[0];
rz(0.42876354) q[1];
sx q[1];
rz(-0.88974297) q[1];
sx q[1];
rz(2.9926328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65229416) q[0];
sx q[0];
rz(-0.51744381) q[0];
sx q[0];
rz(1.6156455) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3690022) q[2];
sx q[2];
rz(-1.861146) q[2];
sx q[2];
rz(-1.5648016) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.91576112) q[1];
sx q[1];
rz(-1.718545) q[1];
sx q[1];
rz(2.3810782) q[1];
rz(0.034986939) q[3];
sx q[3];
rz(-1.9763499) q[3];
sx q[3];
rz(1.4407033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.011977) q[2];
sx q[2];
rz(-2.7836383) q[2];
sx q[2];
rz(0.32148662) q[2];
rz(2.0251515) q[3];
sx q[3];
rz(-0.8588841) q[3];
sx q[3];
rz(1.4139676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1534934) q[0];
sx q[0];
rz(-1.3445925) q[0];
sx q[0];
rz(3.1186812) q[0];
rz(1.5921536) q[1];
sx q[1];
rz(-1.0433082) q[1];
sx q[1];
rz(-1.2124088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3697982) q[0];
sx q[0];
rz(-1.8317458) q[0];
sx q[0];
rz(1.0124932) q[0];
rz(1.171807) q[2];
sx q[2];
rz(-2.3697457) q[2];
sx q[2];
rz(-3.0400624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.089253292) q[1];
sx q[1];
rz(-1.2683378) q[1];
sx q[1];
rz(2.7238728) q[1];
rz(-1.4597662) q[3];
sx q[3];
rz(-2.4250406) q[3];
sx q[3];
rz(2.6727303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93078485) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(-0.99736324) q[2];
rz(-1.8831683) q[3];
sx q[3];
rz(-1.9505898) q[3];
sx q[3];
rz(1.9112126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1579943) q[0];
sx q[0];
rz(-2.4857434) q[0];
sx q[0];
rz(-2.6672145) q[0];
rz(-2.5899218) q[1];
sx q[1];
rz(-0.76849476) q[1];
sx q[1];
rz(-0.65753585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5112572) q[0];
sx q[0];
rz(-1.610582) q[0];
sx q[0];
rz(1.5944832) q[0];
x q[1];
rz(-1.1649407) q[2];
sx q[2];
rz(-1.4144344) q[2];
sx q[2];
rz(0.49977068) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8847981) q[1];
sx q[1];
rz(-1.2024643) q[1];
sx q[1];
rz(0.017466768) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94352888) q[3];
sx q[3];
rz(-2.9637058) q[3];
sx q[3];
rz(-1.2961783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9115596) q[2];
sx q[2];
rz(-0.40731373) q[2];
sx q[2];
rz(-1.3113021) q[2];
rz(-0.71410549) q[3];
sx q[3];
rz(-1.8592535) q[3];
sx q[3];
rz(-0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-3.1353726) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(0.60829341) q[0];
rz(-0.81781203) q[1];
sx q[1];
rz(-1.1759956) q[1];
sx q[1];
rz(2.3086595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.333182) q[0];
sx q[0];
rz(-1.5191374) q[0];
sx q[0];
rz(-0.42148659) q[0];
rz(-pi) q[1];
rz(0.3141381) q[2];
sx q[2];
rz(-1.3262265) q[2];
sx q[2];
rz(0.42942522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3994308) q[1];
sx q[1];
rz(-0.99282904) q[1];
sx q[1];
rz(1.7294243) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5394443) q[3];
sx q[3];
rz(-0.97797457) q[3];
sx q[3];
rz(1.624799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1828764) q[2];
sx q[2];
rz(-0.35280886) q[2];
sx q[2];
rz(-0.58615169) q[2];
rz(-2.2385249) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(0.085845145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387954) q[0];
sx q[0];
rz(-2.0113404) q[0];
sx q[0];
rz(2.197862) q[0];
rz(2.0441652) q[1];
sx q[1];
rz(-2.8891017) q[1];
sx q[1];
rz(-0.15695922) q[1];
rz(1.912276) q[2];
sx q[2];
rz(-1.4394282) q[2];
sx q[2];
rz(-1.9442888) q[2];
rz(0.87370609) q[3];
sx q[3];
rz(-1.9601964) q[3];
sx q[3];
rz(-1.8662069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
