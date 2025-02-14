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
rz(2.4647291) q[0];
sx q[0];
rz(-1.530175) q[0];
sx q[0];
rz(2.8443008) q[0];
rz(-2.3625506) q[1];
sx q[1];
rz(-2.9809597) q[1];
sx q[1];
rz(-1.4359441) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061062) q[0];
sx q[0];
rz(-2.1155042) q[0];
sx q[0];
rz(-2.5404055) q[0];
x q[1];
rz(-1.3117242) q[2];
sx q[2];
rz(-1.4464149) q[2];
sx q[2];
rz(2.4476515) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0983734) q[1];
sx q[1];
rz(-2.3465152) q[1];
sx q[1];
rz(-2.6836392) q[1];
rz(-pi) q[2];
rz(1.3709896) q[3];
sx q[3];
rz(-1.8733896) q[3];
sx q[3];
rz(1.2452919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5982738) q[2];
sx q[2];
rz(-2.7136549) q[2];
sx q[2];
rz(-2.9738026) q[2];
rz(-2.0134036) q[3];
sx q[3];
rz(-1.5690683) q[3];
sx q[3];
rz(-1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0447277) q[0];
sx q[0];
rz(-2.4714578) q[0];
sx q[0];
rz(-2.0665533) q[0];
rz(1.3661522) q[1];
sx q[1];
rz(-1.7551883) q[1];
sx q[1];
rz(-1.8278106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85321097) q[0];
sx q[0];
rz(-0.40069212) q[0];
sx q[0];
rz(-0.59817065) q[0];
rz(-pi) q[1];
rz(2.9552835) q[2];
sx q[2];
rz(-1.7373475) q[2];
sx q[2];
rz(-0.09504091) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3071482) q[1];
sx q[1];
rz(-1.6475369) q[1];
sx q[1];
rz(0.84636023) q[1];
rz(-1.9437378) q[3];
sx q[3];
rz(-2.7881456) q[3];
sx q[3];
rz(1.8693876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65110597) q[2];
sx q[2];
rz(-1.7650812) q[2];
sx q[2];
rz(2.7172078) q[2];
rz(0.73244798) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(1.6781767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3378147) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(2.6388229) q[0];
rz(0.94364199) q[1];
sx q[1];
rz(-2.4991401) q[1];
sx q[1];
rz(-2.3263993) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73683263) q[0];
sx q[0];
rz(-1.5239927) q[0];
sx q[0];
rz(-2.9535049) q[0];
rz(-pi) q[1];
rz(-0.11932245) q[2];
sx q[2];
rz(-1.6762969) q[2];
sx q[2];
rz(1.5169992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6289754) q[1];
sx q[1];
rz(-1.830258) q[1];
sx q[1];
rz(1.8420868) q[1];
rz(2.268154) q[3];
sx q[3];
rz(-2.0202352) q[3];
sx q[3];
rz(-0.12974236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.076685585) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(-2.0150851) q[2];
rz(-1.7846151) q[3];
sx q[3];
rz(-2.6419736) q[3];
sx q[3];
rz(-0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821238) q[0];
sx q[0];
rz(-0.72002763) q[0];
sx q[0];
rz(0.35991392) q[0];
rz(0.24078807) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(2.7684033) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2997595) q[0];
sx q[0];
rz(-1.4275121) q[0];
sx q[0];
rz(-0.10023198) q[0];
rz(-pi) q[1];
rz(2.7794851) q[2];
sx q[2];
rz(-0.59382861) q[2];
sx q[2];
rz(2.1613742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9266921) q[1];
sx q[1];
rz(-0.63758958) q[1];
sx q[1];
rz(-0.59507782) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6945962) q[3];
sx q[3];
rz(-1.742954) q[3];
sx q[3];
rz(-2.3800473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0571664) q[2];
sx q[2];
rz(-1.4320222) q[2];
sx q[2];
rz(-2.2331494) q[2];
rz(1.069979) q[3];
sx q[3];
rz(-1.9570276) q[3];
sx q[3];
rz(2.158304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1835566) q[0];
sx q[0];
rz(-0.89360845) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(2.4519582) q[1];
sx q[1];
rz(-0.90877405) q[1];
sx q[1];
rz(-2.3675809) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2406695) q[0];
sx q[0];
rz(-2.0271447) q[0];
sx q[0];
rz(-0.6610481) q[0];
x q[1];
rz(0.031096259) q[2];
sx q[2];
rz(-0.79582873) q[2];
sx q[2];
rz(-0.21626706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7237295) q[1];
sx q[1];
rz(-2.6681719) q[1];
sx q[1];
rz(-2.0000877) q[1];
x q[2];
rz(-1.6004531) q[3];
sx q[3];
rz(-1.0231442) q[3];
sx q[3];
rz(-0.95049196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8959877) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(0.73949933) q[2];
rz(1.5064404) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(-0.80772775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6226115) q[0];
sx q[0];
rz(-0.7951355) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(-3.0033424) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(1.5672055) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.137321) q[0];
sx q[0];
rz(-1.9372276) q[0];
sx q[0];
rz(1.9204813) q[0];
rz(-pi) q[1];
rz(-0.41344686) q[2];
sx q[2];
rz(-1.8395429) q[2];
sx q[2];
rz(-2.9818802) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3688599) q[1];
sx q[1];
rz(-0.28801685) q[1];
sx q[1];
rz(0.45175987) q[1];
x q[2];
rz(-1.3732143) q[3];
sx q[3];
rz(-1.0302596) q[3];
sx q[3];
rz(-2.5283802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77341998) q[2];
sx q[2];
rz(-1.0680826) q[2];
sx q[2];
rz(2.0424021) q[2];
rz(-2.1558971) q[3];
sx q[3];
rz(-1.516195) q[3];
sx q[3];
rz(-0.43258468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2987591) q[0];
sx q[0];
rz(-2.1465813) q[0];
sx q[0];
rz(0.80879912) q[0];
rz(-0.66566268) q[1];
sx q[1];
rz(-2.1620965) q[1];
sx q[1];
rz(0.98141247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8869276) q[0];
sx q[0];
rz(-1.6279632) q[0];
sx q[0];
rz(-2.1421823) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44156011) q[2];
sx q[2];
rz(-0.54916635) q[2];
sx q[2];
rz(-1.5541981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7715523) q[1];
sx q[1];
rz(-2.0345033) q[1];
sx q[1];
rz(2.9897905) q[1];
rz(-pi) q[2];
rz(-2.2432765) q[3];
sx q[3];
rz(-1.4631464) q[3];
sx q[3];
rz(-0.09569351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.064528) q[2];
sx q[2];
rz(-1.4746102) q[2];
sx q[2];
rz(-2.7725753) q[2];
rz(1.7391694) q[3];
sx q[3];
rz(-0.88715059) q[3];
sx q[3];
rz(1.8203452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8333261) q[0];
sx q[0];
rz(-2.3127191) q[0];
sx q[0];
rz(-2.8114124) q[0];
rz(0.45686832) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(0.59828573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24007455) q[0];
sx q[0];
rz(-0.80722729) q[0];
sx q[0];
rz(0.62348311) q[0];
x q[1];
rz(-3.0706579) q[2];
sx q[2];
rz(-1.1026376) q[2];
sx q[2];
rz(-2.2092961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.48409325) q[1];
sx q[1];
rz(-1.7998905) q[1];
sx q[1];
rz(0.34219663) q[1];
rz(-pi) q[2];
rz(0.65450683) q[3];
sx q[3];
rz(-2.0571567) q[3];
sx q[3];
rz(0.24412046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8699708) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(-3.0769707) q[2];
rz(-1.4501075) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(1.5456642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2794063) q[0];
sx q[0];
rz(-0.33410826) q[0];
sx q[0];
rz(0.52126467) q[0];
rz(2.0817256) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(1.0313787) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17612621) q[0];
sx q[0];
rz(-0.83399189) q[0];
sx q[0];
rz(-1.470572) q[0];
rz(-pi) q[1];
rz(-1.4508444) q[2];
sx q[2];
rz(-1.861683) q[2];
sx q[2];
rz(1.5121258) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.045193521) q[1];
sx q[1];
rz(-2.5303683) q[1];
sx q[1];
rz(-1.0525714) q[1];
rz(-pi) q[2];
rz(-2.6645722) q[3];
sx q[3];
rz(-1.4833772) q[3];
sx q[3];
rz(-1.9893579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3970268) q[2];
sx q[2];
rz(-2.4137745) q[2];
sx q[2];
rz(3.0542206) q[2];
rz(1.9258063) q[3];
sx q[3];
rz(-1.6814503) q[3];
sx q[3];
rz(-2.9714835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755197) q[0];
sx q[0];
rz(-0.89101321) q[0];
sx q[0];
rz(-2.0024306) q[0];
rz(2.6423404) q[1];
sx q[1];
rz(-1.6923994) q[1];
sx q[1];
rz(-1.3334691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3625054) q[0];
sx q[0];
rz(-1.4646069) q[0];
sx q[0];
rz(1.4800735) q[0];
x q[1];
rz(-1.3158445) q[2];
sx q[2];
rz(-1.6189908) q[2];
sx q[2];
rz(-0.15010897) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.002961) q[1];
sx q[1];
rz(-1.8233144) q[1];
sx q[1];
rz(-1.7934402) q[1];
x q[2];
rz(-0.91657775) q[3];
sx q[3];
rz(-2.1567656) q[3];
sx q[3];
rz(2.7295392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.501005) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(2.8583543) q[2];
rz(-1.2011652) q[3];
sx q[3];
rz(-2.4656506) q[3];
sx q[3];
rz(-1.0454319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3129811) q[0];
sx q[0];
rz(-0.77684488) q[0];
sx q[0];
rz(-2.0364398) q[0];
rz(-0.10440566) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(1.3731643) q[2];
sx q[2];
rz(-1.0734954) q[2];
sx q[2];
rz(-1.3263477) q[2];
rz(1.1558044) q[3];
sx q[3];
rz(-0.94528755) q[3];
sx q[3];
rz(-0.62936924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
