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
rz(-2.5928901) q[0];
sx q[0];
rz(-2.2572416) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3724369) q[0];
sx q[0];
rz(-2.7363442) q[0];
sx q[0];
rz(-2.0408819) q[0];
x q[1];
rz(2.5764478) q[2];
sx q[2];
rz(-1.5752715) q[2];
sx q[2];
rz(2.7067513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86815392) q[1];
sx q[1];
rz(-1.3248797) q[1];
sx q[1];
rz(1.6383365) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35688551) q[3];
sx q[3];
rz(-2.2483453) q[3];
sx q[3];
rz(-2.4524636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3551336) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(-0.65594977) q[2];
rz(-1.2077228) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(-0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(0.22878376) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(0.0072335009) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9418966) q[0];
sx q[0];
rz(-1.6640088) q[0];
sx q[0];
rz(1.1774363) q[0];
rz(-2.3184001) q[2];
sx q[2];
rz(-0.41831145) q[2];
sx q[2];
rz(-0.38564607) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.05939535) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(-0.55179623) q[1];
rz(0.90356566) q[3];
sx q[3];
rz(-1.1840608) q[3];
sx q[3];
rz(2.8500593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6146415) q[2];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(-2.8799768) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944377) q[0];
sx q[0];
rz(-1.5804287) q[0];
sx q[0];
rz(-1.5887512) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0171266) q[2];
sx q[2];
rz(-0.83041588) q[2];
sx q[2];
rz(0.24701961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1640472) q[1];
sx q[1];
rz(-0.79037145) q[1];
sx q[1];
rz(2.9708944) q[1];
rz(-pi) q[2];
rz(-1.5393799) q[3];
sx q[3];
rz(-2.0835702) q[3];
sx q[3];
rz(-2.153742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(-2.938081) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(-2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(1.571636) q[0];
rz(-1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.2483695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6682537) q[0];
sx q[0];
rz(-3.0248397) q[0];
sx q[0];
rz(0.97572414) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7522881) q[2];
sx q[2];
rz(-1.776812) q[2];
sx q[2];
rz(2.4455621) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9256546) q[1];
sx q[1];
rz(-2.2207894) q[1];
sx q[1];
rz(-1.4309806) q[1];
rz(-pi) q[2];
rz(-2.4183211) q[3];
sx q[3];
rz(-1.372882) q[3];
sx q[3];
rz(-0.42641446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1127597) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(-0.83703414) q[2];
rz(-1.933243) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(-2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(0.77520448) q[0];
rz(0.40183055) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(-0.90243375) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9426614) q[0];
sx q[0];
rz(-2.5006602) q[0];
sx q[0];
rz(-3.065227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6871704) q[2];
sx q[2];
rz(-1.7917969) q[2];
sx q[2];
rz(-0.99046747) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.24689281) q[1];
sx q[1];
rz(-2.1293853) q[1];
sx q[1];
rz(-0.65319368) q[1];
rz(-pi) q[2];
rz(1.5794472) q[3];
sx q[3];
rz(-0.31673613) q[3];
sx q[3];
rz(-0.99416152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9465785) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(1.9449332) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.3301795) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-0.20176372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70226442) q[0];
sx q[0];
rz(-1.4685255) q[0];
sx q[0];
rz(-3.0754473) q[0];
x q[1];
rz(0.47735881) q[2];
sx q[2];
rz(-1.9743894) q[2];
sx q[2];
rz(-1.8237643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4060496) q[1];
sx q[1];
rz(-0.65945259) q[1];
sx q[1];
rz(2.4156648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1062578) q[3];
sx q[3];
rz(-2.2482276) q[3];
sx q[3];
rz(2.5708452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(2.8743437) q[2];
rz(-2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-0.1593312) q[0];
sx q[0];
rz(-0.057549495) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(-2.1988791) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7563815) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(0.60473196) q[0];
rz(-2.1812181) q[2];
sx q[2];
rz(-2.865961) q[2];
sx q[2];
rz(0.20197091) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47989935) q[1];
sx q[1];
rz(-0.53240314) q[1];
sx q[1];
rz(2.2366621) q[1];
x q[2];
rz(-2.0122635) q[3];
sx q[3];
rz(-0.55674508) q[3];
sx q[3];
rz(0.41551513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(2.7589202) q[2];
rz(1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246178) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-2.8714645) q[0];
rz(-2.5121636) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(-0.28392917) q[1];
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
rz(-pi) q[1];
x q[1];
rz(1.7173041) q[2];
sx q[2];
rz(-2.1657145) q[2];
sx q[2];
rz(-0.3301687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64360196) q[1];
sx q[1];
rz(-1.0712578) q[1];
sx q[1];
rz(-1.040578) q[1];
rz(1.6007401) q[3];
sx q[3];
rz(-0.91120126) q[3];
sx q[3];
rz(-2.3823882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.930687) q[2];
rz(0.25990137) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-0.62817812) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0652086) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(0.22924766) q[0];
rz(-2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.4607666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630336) q[0];
sx q[0];
rz(-1.2872818) q[0];
sx q[0];
rz(-1.589033) q[0];
x q[1];
rz(-1.3912348) q[2];
sx q[2];
rz(-1.944724) q[2];
sx q[2];
rz(-1.1073081) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31635346) q[1];
sx q[1];
rz(-1.5931399) q[1];
sx q[1];
rz(0.52492001) q[1];
rz(0.78482307) q[3];
sx q[3];
rz(-1.8514957) q[3];
sx q[3];
rz(-0.37898889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(-2.4263884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.573695) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(2.8046872) q[0];
rz(-2.9341872) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(2.7609603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56363737) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(1.3689343) q[0];
x q[1];
rz(-2.5194174) q[2];
sx q[2];
rz(-0.84465775) q[2];
sx q[2];
rz(0.13571339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.759728) q[1];
sx q[1];
rz(-1.3761531) q[1];
sx q[1];
rz(-2.0774283) q[1];
rz(2.5467039) q[3];
sx q[3];
rz(-2.4747362) q[3];
sx q[3];
rz(-0.24645933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(-1.1364737) q[2];
rz(0.040955695) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(-1.0271172) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(-0.2480416) q[2];
sx q[2];
rz(-1.7938062) q[2];
sx q[2];
rz(1.2563406) q[2];
rz(-2.9410578) q[3];
sx q[3];
rz(-2.0024828) q[3];
sx q[3];
rz(1.9867292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
