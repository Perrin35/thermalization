OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4857793) q[0];
sx q[0];
rz(-1.9137148) q[0];
sx q[0];
rz(1.9302084) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3762796) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(0.050616654) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.624144) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(1.8271853) q[1];
rz(-pi) q[2];
rz(-2.9458463) q[3];
sx q[3];
rz(-1.9276852) q[3];
sx q[3];
rz(-1.2008047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(1.4663565) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(-1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-2.9557513) q[0];
rz(2.5813685) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-0.21683189) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2730392) q[0];
sx q[0];
rz(-2.2076063) q[0];
sx q[0];
rz(2.7815232) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88044135) q[2];
sx q[2];
rz(-2.464622) q[2];
sx q[2];
rz(1.134269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9358881) q[1];
sx q[1];
rz(-1.3591896) q[1];
sx q[1];
rz(0.88797027) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64036815) q[3];
sx q[3];
rz(-2.159517) q[3];
sx q[3];
rz(-2.0224188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(-2.2089829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55721012) q[0];
sx q[0];
rz(-0.32214468) q[0];
sx q[0];
rz(-1.7412709) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9595397) q[2];
sx q[2];
rz(-1.7943873) q[2];
sx q[2];
rz(1.6857266) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.31308094) q[1];
sx q[1];
rz(-2.5701437) q[1];
sx q[1];
rz(2.9213195) q[1];
rz(-pi) q[2];
rz(-2.1266537) q[3];
sx q[3];
rz(-1.9564637) q[3];
sx q[3];
rz(2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(1.0926584) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(-0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(-1.6436228) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.362975) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(2.1182563) q[0];
rz(1.3802714) q[2];
sx q[2];
rz(-1.6263783) q[2];
sx q[2];
rz(-1.3560825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8343463) q[1];
sx q[1];
rz(-1.3386968) q[1];
sx q[1];
rz(-2.8744065) q[1];
rz(0.46521503) q[3];
sx q[3];
rz(-1.1606996) q[3];
sx q[3];
rz(2.0801534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3952289) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(-0.70181075) q[2];
rz(2.3102405) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(-2.3262614) q[0];
rz(-1.6197846) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-1.048208) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3807555) q[0];
sx q[0];
rz(-1.8341944) q[0];
sx q[0];
rz(-2.8834228) q[0];
x q[1];
rz(2.4900715) q[2];
sx q[2];
rz(-1.5805575) q[2];
sx q[2];
rz(2.5196645) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.293562) q[1];
sx q[1];
rz(-1.3333496) q[1];
sx q[1];
rz(-2.8192987) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41042495) q[3];
sx q[3];
rz(-0.99567185) q[3];
sx q[3];
rz(2.9068973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(-0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734171) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-3.0117603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8311365) q[0];
sx q[0];
rz(-1.2015011) q[0];
sx q[0];
rz(0.62311689) q[0];
rz(-0.36826276) q[2];
sx q[2];
rz(-1.0266745) q[2];
sx q[2];
rz(0.95168176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.063720238) q[1];
sx q[1];
rz(-1.93396) q[1];
sx q[1];
rz(-2.5329464) q[1];
rz(-pi) q[2];
rz(0.54883212) q[3];
sx q[3];
rz(-1.4026814) q[3];
sx q[3];
rz(-2.7041534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(2.9373346) q[2];
rz(-1.2060818) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(-2.4553305) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26353729) q[0];
sx q[0];
rz(-1.1362846) q[0];
sx q[0];
rz(-2.6146019) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1967728) q[2];
sx q[2];
rz(-0.84569028) q[2];
sx q[2];
rz(-0.041989728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2974907) q[1];
sx q[1];
rz(-2.9417848) q[1];
sx q[1];
rz(1.954345) q[1];
rz(-pi) q[2];
rz(-2.5043082) q[3];
sx q[3];
rz(-2.6316959) q[3];
sx q[3];
rz(-2.2989458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(3.1398204) q[2];
rz(-2.5799675) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(-1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(1.6954533) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(1.5015645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1743463) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(-1.0780225) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1649706) q[2];
sx q[2];
rz(-0.067194447) q[2];
sx q[2];
rz(-0.42025987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20128076) q[1];
sx q[1];
rz(-2.8021325) q[1];
sx q[1];
rz(-2.8708354) q[1];
rz(-pi) q[2];
rz(-2.4092259) q[3];
sx q[3];
rz(-0.77419705) q[3];
sx q[3];
rz(-2.458651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(-1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(2.8222728) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.2040899) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8512745) q[0];
sx q[0];
rz(-1.9418678) q[0];
sx q[0];
rz(2.0593658) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7792764) q[2];
sx q[2];
rz(-3*pi/13) q[2];
sx q[2];
rz(2.97646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0723567) q[1];
sx q[1];
rz(-1.623739) q[1];
sx q[1];
rz(-2.3318021) q[1];
rz(-pi) q[2];
rz(-1.5185235) q[3];
sx q[3];
rz(-1.4613323) q[3];
sx q[3];
rz(-1.048552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7982771) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(-0.19432755) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1949085) q[0];
sx q[0];
rz(-2.2241728) q[0];
sx q[0];
rz(0.16434591) q[0];
rz(0.98722234) q[2];
sx q[2];
rz(-0.7910896) q[2];
sx q[2];
rz(2.1949878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.918805) q[1];
sx q[1];
rz(-0.58090392) q[1];
sx q[1];
rz(-1.3851628) q[1];
rz(-pi) q[2];
rz(-2.4284027) q[3];
sx q[3];
rz(-1.4048368) q[3];
sx q[3];
rz(0.99456577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.4355961) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(-0.96799093) q[2];
sx q[2];
rz(-1.597076) q[2];
sx q[2];
rz(1.1521641) q[2];
rz(0.070449645) q[3];
sx q[3];
rz(-1.9000713) q[3];
sx q[3];
rz(-2.6303359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
