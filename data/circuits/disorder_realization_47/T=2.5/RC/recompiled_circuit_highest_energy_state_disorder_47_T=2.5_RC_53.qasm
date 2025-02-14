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
rz(0.53606755) q[0];
sx q[0];
rz(7.4458691) q[0];
sx q[0];
rz(9.8557628) q[0];
rz(-2.7886136) q[1];
sx q[1];
rz(-1.2231491) q[1];
sx q[1];
rz(-0.42981848) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0551222) q[0];
sx q[0];
rz(-1.6608547) q[0];
sx q[0];
rz(-0.49205972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.016021803) q[2];
sx q[2];
rz(-1.1384769) q[2];
sx q[2];
rz(-0.95903083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42001549) q[1];
sx q[1];
rz(-2.7215981) q[1];
sx q[1];
rz(0.53424044) q[1];
rz(2.6453703) q[3];
sx q[3];
rz(-1.2321951) q[3];
sx q[3];
rz(2.7994775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.063804403) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(-1.0112313) q[2];
rz(1.0412591) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(2.7019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1083199) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(-2.394115) q[0];
rz(-0.29391089) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(-2.8864313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9440112) q[0];
sx q[0];
rz(-1.4048049) q[0];
sx q[0];
rz(0.98910467) q[0];
rz(-pi) q[1];
rz(2.8361144) q[2];
sx q[2];
rz(-1.3925838) q[2];
sx q[2];
rz(-0.62159789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9772759) q[1];
sx q[1];
rz(-2.212388) q[1];
sx q[1];
rz(2.5733828) q[1];
x q[2];
rz(-1.3142845) q[3];
sx q[3];
rz(-0.64033901) q[3];
sx q[3];
rz(-0.72968715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1110288) q[2];
sx q[2];
rz(-0.054628987) q[2];
sx q[2];
rz(0.89152208) q[2];
rz(-0.24957481) q[3];
sx q[3];
rz(-2.2587743) q[3];
sx q[3];
rz(0.3961302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7340649) q[0];
sx q[0];
rz(-1.1363109) q[0];
sx q[0];
rz(-3.0974645) q[0];
rz(1.2491501) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(-2.6557907) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51903782) q[0];
sx q[0];
rz(-1.7960894) q[0];
sx q[0];
rz(-1.5908498) q[0];
rz(-0.13834149) q[2];
sx q[2];
rz(-1.7560343) q[2];
sx q[2];
rz(-1.1323765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1959744) q[1];
sx q[1];
rz(-0.83688078) q[1];
sx q[1];
rz(-2.3952146) q[1];
x q[2];
rz(2.0030735) q[3];
sx q[3];
rz(-1.2069993) q[3];
sx q[3];
rz(-1.0897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8140063) q[2];
sx q[2];
rz(-2.0531211) q[2];
sx q[2];
rz(2.4468454) q[2];
rz(-0.57573777) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(1.7339138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1908252) q[0];
sx q[0];
rz(-2.8468843) q[0];
sx q[0];
rz(-2.8488391) q[0];
rz(-0.5689019) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(0.84691602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900565) q[0];
sx q[0];
rz(-2.7256294) q[0];
sx q[0];
rz(-1.7638168) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0693502) q[2];
sx q[2];
rz(-2.2390167) q[2];
sx q[2];
rz(-2.7492858) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21601115) q[1];
sx q[1];
rz(-2.8329509) q[1];
sx q[1];
rz(0.68283783) q[1];
rz(-pi) q[2];
x q[2];
rz(0.067518721) q[3];
sx q[3];
rz(-1.3775) q[3];
sx q[3];
rz(-0.75261469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7920821) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(2.9328031) q[2];
rz(-0.71792349) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(-0.94055241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1013041) q[0];
sx q[0];
rz(-2.6191819) q[0];
sx q[0];
rz(-2.8029602) q[0];
rz(2.4385117) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(-0.83831659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378483) q[0];
sx q[0];
rz(-0.5068584) q[0];
sx q[0];
rz(0.54308191) q[0];
x q[1];
rz(-1.3569068) q[2];
sx q[2];
rz(-0.75816064) q[2];
sx q[2];
rz(-1.2616829) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.053239659) q[1];
sx q[1];
rz(-0.87118282) q[1];
sx q[1];
rz(-3.0979473) q[1];
rz(-pi) q[2];
rz(-0.60689599) q[3];
sx q[3];
rz(-1.9268774) q[3];
sx q[3];
rz(2.0739561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4579953) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(-1.0028769) q[2];
rz(2.9966127) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(0.063044757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.018464) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(2.1783094) q[0];
rz(0.18114289) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(1.9021665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9267777) q[0];
sx q[0];
rz(-0.96450761) q[0];
sx q[0];
rz(-1.1112332) q[0];
rz(-pi) q[1];
rz(2.6341637) q[2];
sx q[2];
rz(-2.6118738) q[2];
sx q[2];
rz(1.5629753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9628004) q[1];
sx q[1];
rz(-1.9427248) q[1];
sx q[1];
rz(0.40403251) q[1];
x q[2];
rz(-1.7224947) q[3];
sx q[3];
rz(-1.6983508) q[3];
sx q[3];
rz(-0.90326819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.88508254) q[2];
sx q[2];
rz(-0.91826597) q[2];
sx q[2];
rz(2.2086842) q[2];
rz(1.4126623) q[3];
sx q[3];
rz(-1.7191073) q[3];
sx q[3];
rz(-2.9322114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949718) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(-0.80286655) q[0];
rz(0.74202263) q[1];
sx q[1];
rz(-1.6525533) q[1];
sx q[1];
rz(2.7313357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0502346) q[0];
sx q[0];
rz(-1.0438393) q[0];
sx q[0];
rz(-1.1731253) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1226317) q[2];
sx q[2];
rz(-0.76803128) q[2];
sx q[2];
rz(2.1036069) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2563808) q[1];
sx q[1];
rz(-1.4090013) q[1];
sx q[1];
rz(0.075116861) q[1];
rz(1.0954082) q[3];
sx q[3];
rz(-2.1719526) q[3];
sx q[3];
rz(0.22829311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86847574) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(-0.22129076) q[2];
rz(-0.41915974) q[3];
sx q[3];
rz(-2.6594888) q[3];
sx q[3];
rz(2.6656849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047423) q[0];
sx q[0];
rz(-0.99822799) q[0];
sx q[0];
rz(1.6118443) q[0];
rz(1.8086241) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(1.6882247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71163342) q[0];
sx q[0];
rz(-2.2471513) q[0];
sx q[0];
rz(2.3948689) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55744268) q[2];
sx q[2];
rz(-0.22391437) q[2];
sx q[2];
rz(-0.73665184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9295974) q[1];
sx q[1];
rz(-0.72337615) q[1];
sx q[1];
rz(0.22774793) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2923041) q[3];
sx q[3];
rz(-1.7207816) q[3];
sx q[3];
rz(1.5962792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2754485) q[2];
sx q[2];
rz(-0.25442213) q[2];
sx q[2];
rz(3.1308543) q[2];
rz(2.8209316) q[3];
sx q[3];
rz(-1.843822) q[3];
sx q[3];
rz(2.5519154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7084932) q[0];
sx q[0];
rz(-2.1864317) q[0];
sx q[0];
rz(0.47472111) q[0];
rz(-2.8482598) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(1.6620103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5005762) q[0];
sx q[0];
rz(-1.5307284) q[0];
sx q[0];
rz(3.1047719) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8117254) q[2];
sx q[2];
rz(-2.2325667) q[2];
sx q[2];
rz(2.2994201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7254527) q[1];
sx q[1];
rz(-1.801277) q[1];
sx q[1];
rz(-2.6621755) q[1];
rz(2.1988499) q[3];
sx q[3];
rz(-0.30170479) q[3];
sx q[3];
rz(0.42886558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.097229615) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(-2.7105159) q[2];
rz(-2.9512682) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(-0.88461191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97656074) q[0];
sx q[0];
rz(-0.29104069) q[0];
sx q[0];
rz(-0.2188368) q[0];
rz(0.95895514) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(-1.3921907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4138654) q[0];
sx q[0];
rz(-1.6221592) q[0];
sx q[0];
rz(-0.19330175) q[0];
rz(1.7151681) q[2];
sx q[2];
rz(-1.1544268) q[2];
sx q[2];
rz(-1.8954111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4569623) q[1];
sx q[1];
rz(-0.70856386) q[1];
sx q[1];
rz(2.674391) q[1];
rz(-0.95646071) q[3];
sx q[3];
rz(-1.9298565) q[3];
sx q[3];
rz(-0.39310716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86768156) q[2];
sx q[2];
rz(-0.76649222) q[2];
sx q[2];
rz(-1.8217746) q[2];
rz(-2.6093318) q[3];
sx q[3];
rz(-1.3479193) q[3];
sx q[3];
rz(-1.5490612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821447) q[0];
sx q[0];
rz(-1.5503217) q[0];
sx q[0];
rz(0.4009608) q[0];
rz(-2.9136912) q[1];
sx q[1];
rz(-2.0961998) q[1];
sx q[1];
rz(-1.6963522) q[1];
rz(-2.5980742) q[2];
sx q[2];
rz(-0.47558161) q[2];
sx q[2];
rz(-2.4408792) q[2];
rz(2.1501272) q[3];
sx q[3];
rz(-1.7064863) q[3];
sx q[3];
rz(1.9319921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
