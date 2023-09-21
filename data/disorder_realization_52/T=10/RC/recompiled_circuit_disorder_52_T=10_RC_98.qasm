OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(6.8318879) q[0];
sx q[0];
rz(5.3988342) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(1.6391099) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76915574) q[0];
sx q[0];
rz(-0.40524846) q[0];
sx q[0];
rz(-1.1007108) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0083562334) q[2];
sx q[2];
rz(-2.5764321) q[2];
sx q[2];
rz(-1.9985808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2734387) q[1];
sx q[1];
rz(-1.3248797) q[1];
sx q[1];
rz(1.5032561) q[1];
x q[2];
rz(0.35688551) q[3];
sx q[3];
rz(-2.2483453) q[3];
sx q[3];
rz(-0.68912904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78645906) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(2.4856429) q[2];
rz(-1.9338699) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2475964) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(-0.43352747) q[0];
rz(0.22878376) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(3.1343592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4097071) q[0];
sx q[0];
rz(-1.1792372) q[0];
sx q[0];
rz(-0.10086985) q[0];
rz(1.2556633) q[2];
sx q[2];
rz(-1.2909781) q[2];
sx q[2];
rz(1.8880106) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65456395) q[1];
sx q[1];
rz(-1.9332758) q[1];
sx q[1];
rz(-1.3351721) q[1];
rz(0.47827999) q[3];
sx q[3];
rz(-0.96049958) q[3];
sx q[3];
rz(1.5680716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5269512) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(-2.4439404) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(-0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4678629) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(0.26161584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3221489) q[0];
sx q[0];
rz(-1.5804287) q[0];
sx q[0];
rz(1.5528414) q[0];
rz(-pi) q[1];
rz(-3.0171266) q[2];
sx q[2];
rz(-0.83041588) q[2];
sx q[2];
rz(2.894573) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6690327) q[1];
sx q[1];
rz(-1.6918039) q[1];
sx q[1];
rz(-2.3585412) q[1];
rz(-pi) q[2];
rz(-2.6286078) q[3];
sx q[3];
rz(-1.5981711) q[3];
sx q[3];
rz(-2.5740636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(-2.938081) q[2];
rz(-0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(-2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1638284) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.571636) q[0];
rz(2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.2483695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0715863) q[0];
sx q[0];
rz(-1.6674111) q[0];
sx q[0];
rz(-3.0759401) q[0];
rz(-2.4291971) q[2];
sx q[2];
rz(-2.8678896) q[2];
sx q[2];
rz(1.7143539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9256546) q[1];
sx q[1];
rz(-0.92080322) q[1];
sx q[1];
rz(1.4309806) q[1];
rz(-pi) q[2];
rz(-2.4183211) q[3];
sx q[3];
rz(-1.7687106) q[3];
sx q[3];
rz(0.42641446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0288329) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(1.933243) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0060624881) q[0];
sx q[0];
rz(-1.0536138) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(0.90243375) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1989312) q[0];
sx q[0];
rz(-2.5006602) q[0];
sx q[0];
rz(-3.065227) q[0];
rz(-2.6685153) q[2];
sx q[2];
rz(-2.6396751) q[2];
sx q[2];
rz(-2.9830473) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.21197) q[1];
sx q[1];
rz(-2.3096482) q[1];
sx q[1];
rz(2.3421939) q[1];
rz(-pi) q[2];
rz(0.0028354672) q[3];
sx q[3];
rz(-1.2540725) q[3];
sx q[3];
rz(-2.1383274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9465785) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(1.6992016) q[3];
sx q[3];
rz(-1.8701575) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(-2.6384171) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(0.20176372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4393282) q[0];
sx q[0];
rz(-1.4685255) q[0];
sx q[0];
rz(3.0754473) q[0];
rz(-2.0189507) q[2];
sx q[2];
rz(-1.1345703) q[2];
sx q[2];
rz(-3.089038) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4060496) q[1];
sx q[1];
rz(-2.4821401) q[1];
sx q[1];
rz(2.4156648) q[1];
x q[2];
rz(1.1062578) q[3];
sx q[3];
rz(-0.89336508) q[3];
sx q[3];
rz(-2.5708452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(-0.26724896) q[2];
rz(2.3184508) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(-2.2657623) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(0.057549495) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(-2.1988791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71390401) q[0];
sx q[0];
rz(-1.1770935) q[0];
sx q[0];
rz(-2.4967771) q[0];
rz(-pi) q[1];
rz(2.9808688) q[2];
sx q[2];
rz(-1.3458999) q[2];
sx q[2];
rz(-0.42663923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2191636) q[1];
sx q[1];
rz(-1.1601829) q[1];
sx q[1];
rz(-0.34904042) q[1];
rz(-pi) q[2];
x q[2];
rz(1.058217) q[3];
sx q[3];
rz(-1.7985385) q[3];
sx q[3];
rz(-2.3678569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1192347) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(0.38267246) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7169749) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-2.8714645) q[0];
rz(2.5121636) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(2.8576635) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086872452) q[0];
sx q[0];
rz(-0.76612872) q[0];
sx q[0];
rz(0.92673577) q[0];
rz(-0.21247272) q[2];
sx q[2];
rz(-0.61056911) q[2];
sx q[2];
rz(0.072710466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8995081) q[1];
sx q[1];
rz(-2.429932) q[1];
sx q[1];
rz(-2.3942024) q[1];
x q[2];
rz(-1.6007401) q[3];
sx q[3];
rz(-0.91120126) q[3];
sx q[3];
rz(-0.75920446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5876864) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(1.930687) q[2];
rz(-0.25990137) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07638409) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(2.912345) q[0];
rz(-2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.4607666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630336) q[0];
sx q[0];
rz(-1.8543108) q[0];
sx q[0];
rz(-1.589033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7503579) q[2];
sx q[2];
rz(-1.1968687) q[2];
sx q[2];
rz(2.0342846) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2673805) q[1];
sx q[1];
rz(-1.0460209) q[1];
sx q[1];
rz(1.5966148) q[1];
x q[2];
rz(-0.38735729) q[3];
sx q[3];
rz(-2.318317) q[3];
sx q[3];
rz(2.2203317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71172697) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(-2.7097278) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(-2.4263884) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.573695) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(0.33690548) q[0];
rz(2.9341872) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(2.7609603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56363737) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(1.7726583) q[0];
rz(-0.74110465) q[2];
sx q[2];
rz(-2.0217102) q[2];
sx q[2];
rz(1.2620743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14643529) q[1];
sx q[1];
rz(-2.6019115) q[1];
sx q[1];
rz(1.1848918) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5638644) q[3];
sx q[3];
rz(-1.9247705) q[3];
sx q[3];
rz(2.3059394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(-2.005119) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.3142746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(1.0271172) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(0.2480416) q[2];
sx q[2];
rz(-1.3477865) q[2];
sx q[2];
rz(-1.8852521) q[2];
rz(-0.20053486) q[3];
sx q[3];
rz(-1.1391098) q[3];
sx q[3];
rz(-1.1548635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
